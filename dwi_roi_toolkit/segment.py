#!/usr/bin/env python3
"""
肺占位 ROI 半自动全体积勾画 (DWI 高信号病灶)

设计前提: 靶区由医生指定种子点, 脚本负责把种子长成全体积 ROI。
脚本不做"哪个是病灶"的判断 —— 那是阅片医生的事。

用法:
    # 医生在 ITK-SNAP 里读出病灶中心体素坐标 (x y z), 可给多个种子
    python segment_roi.py --image DWI_b800.nii.gz --seed 210 205 44 --out out_dir
    # 已有粗 mask 时可用 --init-mask 代替种子
    python segment_roi.py --image DWI_b800.nii.gz --init-mask rough.nii.gz --out out_dir
    # 有 ADC 图时一并统计
    python segment_roi.py --image DWI_b800.nii.gz --seed 210 205 44 --adc ADC.nii.gz --out out_dir

流程:
    各向异性扩散去噪 -> ConfidenceConnected 区域生长 -> 形状约束 level set 细化
    -> 形态学补洞去毛刺 -> 取种子所在连通域 -> 输出 label + 统计表
"""
import os, json, argparse
import numpy as np
import SimpleITK as sitk


def denoise(img):
    """曲率各向异性扩散: DWI 信噪比低, 直接生长会漏出去; 该滤波保边去噪"""
    f = sitk.CurvatureAnisotropicDiffusionImageFilter()
    f.SetTimeStep(0.0625)       # 3D 稳定性上限是 1/2^(dim+1)
    f.SetNumberOfIterations(5)
    f.SetConductanceParameter(2.0)
    return f.Execute(sitk.Cast(img, sitk.sitkFloat32))


def grow(img, seeds, multiplier, radius, iterations):
    """ConfidenceConnected: 用种子邻域的均值±k*标准差自适应定阈值, 比固定阈值稳"""
    f = sitk.ConfidenceConnectedImageFilter()
    f.SetSeedList([tuple(int(v) for v in s) for s in seeds])
    f.SetMultiplier(multiplier)          # k, 越大长得越开
    f.SetInitialNeighborhoodRadius(radius)
    f.SetNumberOfIterations(iterations)  # 每轮用当前区域重估均值方差
    f.SetReplaceValue(1)
    return f.Execute(img)


def refine_levelset(img, init_mask, curvature=1.0, propagation=1.0, n_iter=200):
    """
    阈值 level set 细化边界: 区域生长的边界锯齿明显, level set 用曲率项磨平,
    同时靠强度阈值把边界钉在病灶-肺实质交界上。
    """
    arr = sitk.GetArrayFromImage(img)
    m = sitk.GetArrayFromImage(init_mask).astype(bool)
    if m.sum() < 10:
        return init_mask
    vals = arr[m]
    lo, hi = float(np.percentile(vals, 5)), float(arr.max())

    # level set 要求初始为有符号距离场, 内部为负
    init = sitk.SignedMaurerDistanceMap(init_mask, insideIsPositive=False,
                                        squaredDistance=False, useImageSpacing=True)
    ls = sitk.ThresholdSegmentationLevelSetImageFilter()
    ls.SetLowerThreshold(lo)
    ls.SetUpperThreshold(hi)
    ls.SetMaximumRMSError(0.02)
    ls.SetNumberOfIterations(n_iter)
    ls.SetCurvatureScaling(curvature)        # 平滑力
    ls.SetPropagationScaling(propagation)    # 膨胀力
    out = ls.Execute(sitk.Cast(init, sitk.sitkFloat32), sitk.Cast(img, sitk.sitkFloat32))
    refined = sitk.Cast(out < 0, sitk.sitkUInt8)

    # 安全阀: 小病灶上曲率项可能把轮廓整个收没, 或者阈值太松导致爆炸性外扩。
    # 两种情况都退回区域生长结果, 宁可边界毛糙也不能交一个空的或糊掉的 ROI。
    n_raw = int(sitk.GetArrayFromImage(init_mask).sum())
    n_new = int(sitk.GetArrayFromImage(refined).sum())
    if n_new < 0.2 * n_raw or n_new > 5.0 * n_raw:
        import warnings
        warnings.warn(
            f"level set 结果异常 ({n_raw} -> {n_new} 体素), 已退回区域生长结果。"
            f"可调 --curvature / --propagation, 或加 --no-levelset。",
            RuntimeWarning, stacklevel=2)
        return init_mask
    return refined


def cleanup(mask, seeds, min_voxels=10):
    """补洞 + 开运算去毛刺 + 只保留种子所在的那个连通域"""
    mask = sitk.BinaryFillhole(mask)
    mask = sitk.BinaryMorphologicalOpening(mask, [1, 1, 1], sitk.sitkBall)
    cc = sitk.RelabelComponent(sitk.ConnectedComponent(mask), minimumObjectSize=min_voxels)
    arr = sitk.GetArrayFromImage(cc)
    keep = set()
    for x, y, z in seeds:
        # SimpleITK 索引是 (x,y,z), numpy 数组是 [z,y,x]
        lbl = int(arr[int(z), int(y), int(x)])
        if lbl > 0:
            keep.add(lbl)
    if not keep:                       # 种子落在了被清掉的小碎块上, 退回最大连通域
        keep = {1} if arr.max() >= 1 else set()
    out = np.isin(arr, list(keep)).astype(np.uint8)
    res = sitk.GetImageFromArray(out)
    res.CopyInformation(mask)          # 关键: 保住 spacing/origin/direction
    return res


def stats(image, mask, adc=None):
    """体积 + 信号统计; 有 ADC 时给 ADCmean/min, 前瞻研究常用这几个指标"""
    f = sitk.LabelStatisticsImageFilter()
    f.Execute(sitk.Cast(image, sitk.sitkFloat32), sitk.Cast(mask, sitk.sitkUInt8))
    if 1 not in f.GetLabels():
        return {"error": "ROI 为空, 请调整 --multiplier 或换种子点"}
    vox_mm3 = float(np.prod(mask.GetSpacing()))
    n = int(f.GetCount(1))
    out = {
        "voxels": n,
        "voxel_volume_mm3": round(vox_mm3, 4),
        "volume_mm3": round(n * vox_mm3, 2),
        "volume_ml": round(n * vox_mm3 / 1000.0, 4),
        "signal_mean": round(float(f.GetMean(1)), 3),
        "signal_sd": round(float(f.GetSigma(1)), 3),
        "signal_min": round(float(f.GetMinimum(1)), 3),
        "signal_max": round(float(f.GetMaximum(1)), 3),
    }
    if adc is not None:
        a = sitk.LabelStatisticsImageFilter()
        adc_r = sitk.Resample(adc, mask, sitk.Transform(), sitk.sitkLinear, 0.0, sitk.sitkFloat32)
        a.Execute(adc_r, sitk.Cast(mask, sitk.sitkUInt8))
        if 1 in a.GetLabels():
            out.update({
                "ADC_mean": round(float(a.GetMean(1)), 3),
                "ADC_sd": round(float(a.GetSigma(1)), 3),
                "ADC_min": round(float(a.GetMinimum(1)), 3),
            })
    return out


def build_parser(ap):
    ap.add_argument("--image", required=True, help="分割所用图像 (通常是高 b 值 DWI)")
    ap.add_argument("--seed", nargs=3, type=int, action="append",
                    help="种子体素坐标 x y z (ITK-SNAP 光标处的 index), 可重复给多个")
    ap.add_argument("--init-mask", help="已有粗 mask, 代替种子点")
    ap.add_argument("--adc", help="ADC 图, 用于统计")
    ap.add_argument("--out", required=True)
    ap.add_argument("--multiplier", type=float, default=2.5, help="生长松紧度, 漏出去就调小")
    ap.add_argument("--radius", type=int, default=2)
    ap.add_argument("--iterations", type=int, default=5)
    ap.add_argument("--no-levelset", action="store_true")
    ap.add_argument("--propagation", type=float, default=1.0,
                    help="level set 膨胀力, 边界外扩就调小")
    ap.add_argument("--curvature", type=float, default=1.0,
                    help="level set 平滑力, 小病灶被磨没就调小")
    return ap


def run(a):
    if not a.seed and not a.init_mask:
        raise SystemExit("必须给 --seed 或 --init-mask")
    os.makedirs(a.out, exist_ok=True)

    img = sitk.ReadImage(a.image)
    sm = denoise(img)

    if a.init_mask:
        raw = sitk.Cast(sitk.ReadImage(a.init_mask) > 0, sitk.sitkUInt8)
        arr = sitk.GetArrayFromImage(raw)
        zs, ys, xs = np.where(arr > 0)
        seeds = [(int(xs.mean()), int(ys.mean()), int(zs.mean()))]  # 质心当种子用于保连通域
    else:
        seeds = a.seed
        raw = grow(sm, seeds, a.multiplier, a.radius, a.iterations)

    mask = raw if a.no_levelset else refine_levelset(
        sm, raw, curvature=a.curvature, propagation=a.propagation)
    mask = cleanup(mask, seeds)

    # 输出 uint8 label —— ITK-SNAP 的 Segmentation 要求整数类型
    label_path = os.path.join(a.out, "label.nii.gz")
    sitk.WriteImage(sitk.Cast(mask, sitk.sitkUInt8), label_path, True)
    sitk.WriteImage(img, os.path.join(a.out, "image.nii.gz"), True)

    s = stats(img, mask, sitk.ReadImage(a.adc, sitk.sitkFloat32) if a.adc else None)
    s["seeds"] = seeds
    s["source_image"] = os.path.basename(a.image)
    with open(os.path.join(a.out, "roi_stats.json"), "w", encoding="utf-8") as f:
        json.dump(s, f, ensure_ascii=False, indent=2)

    print(json.dumps(s, ensure_ascii=False, indent=2))
    print(f"\nlabel -> {label_path}")
    print("ITK-SNAP: File > Open Main Image 选 image.nii.gz, 再 Segmentation > Open Segmentation 选 label.nii.gz")


if __name__ == "__main__":
    run(build_parser(argparse.ArgumentParser()).parse_args())
