#!/usr/bin/env python3
"""
非刚性配准: 把 DWI (moving) 配到结构像 (fixed, 如 T2/T1/CT)

用法:
    python register.py --fixed T2.nii.gz --moving DWI_b800.nii.gz --out out_dir \
                       [--also-warp ADC.nii.gz label.nii.gz] [--mode bspline|demons]

策略 (胸部呼吸运动形变大, 分三级由粗到细):
    1. 质心初始化
    2. 刚性 (Versor3D) —— 纠正整体平移旋转
    3. B-spline 自由形变 —— 纠正呼吸/心跳导致的局部形变
    多分辨率金字塔 + 高斯平滑, 避免落进局部极小值。

输出:
    moved.nii.gz        配准后的 DWI
    transform.tfm       复合变换 (可复用到 ADC / label / 其他 b 值)
    displacement.nii.gz 形变场 (供质控)
    checker.nii.gz      棋盘格对比图 (ITK-SNAP 里看对齐效果)
"""
import os, argparse
import SimpleITK as sitk


def log(reg):
    print(f"  level iter={reg.GetOptimizerIteration():3d}  metric={reg.GetMetricValue():.5f}")


def make_rigid(fixed, moving):
    """第一步: 刚性配准, 结果作为形变配准的初值"""
    init = sitk.CenteredTransformInitializer(
        fixed, moving, sitk.VersorRigid3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY)

    reg = sitk.ImageRegistrationMethod()
    # 互信息: DWI 与 T2/CT 对比度不同, 属多模态, 不能用 MSE
    reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.20, seed=42)   # 固定 seed 保证可复现
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsRegularStepGradientDescent(
        learningRate=1.0, minStep=1e-4, numberOfIterations=200,
        gradientMagnitudeTolerance=1e-6)
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel([4, 2, 1])
    reg.SetSmoothingSigmasPerLevel([2, 1, 0])
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    reg.SetInitialTransform(init, inPlace=False)
    reg.AddCommand(sitk.sitkIterationEvent, lambda: log(reg))

    print("[1/2] 刚性配准 ...")
    return reg.Execute(fixed, moving)


def make_bspline(fixed, moving, initial):
    """第二步: B-spline 自由形变"""
    # 控制点网格 ~ 每 30mm 一个, 胸部这个尺度能吃住呼吸形变又不至于过拟合
    mesh = [max(2, int(sz * sp / 30.0)) for sz, sp in zip(fixed.GetSize(), fixed.GetSpacing())]
    bspline = sitk.BSplineTransformInitializer(fixed, mesh, order=3)

    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.20, seed=42)
    reg.SetInterpolator(sitk.sitkLinear)
    # LBFGSB 对 B-spline 系数这种高维参数收敛稳定
    reg.SetOptimizerAsLBFGSB(gradientConvergenceTolerance=1e-5,
                             numberOfIterations=150, maximumNumberOfCorrections=5)
    reg.SetShrinkFactorsPerLevel([4, 2, 1])
    reg.SetSmoothingSigmasPerLevel([2, 1, 0])
    reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    reg.SetMovingInitialTransform(initial)          # 刚性结果作为初值
    reg.SetInitialTransformAsBSpline(bspline, inPlace=True, scaleFactors=[1, 2, 4])
    reg.AddCommand(sitk.sitkIterationEvent, lambda: log(reg))

    print("[2/2] B-spline 形变配准 ...")
    out = reg.Execute(fixed, moving)
    composite = sitk.CompositeTransform([initial, out])
    return composite


def make_demons(fixed, moving, initial):
    """备选: Diffeomorphic Demons, 单模态(如 DWI 不同 b 值之间)更快更稳"""
    pre = sitk.Resample(moving, fixed, initial, sitk.sitkLinear, 0.0, moving.GetPixelID())
    # 直方图匹配是 Demons 的前提, 它假设同一组织灰度相同
    pre = sitk.HistogramMatching(sitk.Cast(pre, sitk.sitkFloat32),
                                 sitk.Cast(fixed, sitk.sitkFloat32),
                                 numberOfHistogramLevels=1024, numberOfMatchPoints=7)
    demons = sitk.FastSymmetricForcesDemonsRegistrationFilter()
    demons.SetNumberOfIterations(100)
    demons.SetStandardDeviations(1.5)
    field = demons.Execute(sitk.Cast(fixed, sitk.sitkFloat32), pre)
    return sitk.CompositeTransform([initial, sitk.DisplacementFieldTransform(field)])


def build_parser(ap):
    ap.add_argument("--fixed", required=True, help="参考图 (结构像)")
    ap.add_argument("--moving", required=True, help="待配准图 (DWI)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", default="bspline", choices=["bspline", "demons"])
    ap.add_argument("--also-warp", nargs="*", default=[],
                    help="用同一变换重采样的其他图 (ADC / 其他 b 值 / 已有 label)")
    return ap


def run(a):
    os.makedirs(a.out, exist_ok=True)

    fixed = sitk.Cast(sitk.ReadImage(a.fixed), sitk.sitkFloat32)
    moving = sitk.Cast(sitk.ReadImage(a.moving), sitk.sitkFloat32)

    rigid = make_rigid(fixed, moving)
    tx = make_bspline(fixed, moving, rigid) if a.mode == "bspline" else make_demons(fixed, moving, rigid)

    sitk.WriteTransform(tx, os.path.join(a.out, "transform.tfm"))
    moved = sitk.Resample(moving, fixed, tx, sitk.sitkBSpline, 0.0, sitk.sitkFloat32)
    sitk.WriteImage(moved, os.path.join(a.out, "moved.nii.gz"), True)

    # 形变场 + 棋盘格, 交付时给对方做目视质控
    field = sitk.TransformToDisplacementField(
        tx, sitk.sitkVectorFloat64, fixed.GetSize(), fixed.GetOrigin(),
        fixed.GetSpacing(), fixed.GetDirection())
    sitk.WriteImage(field, os.path.join(a.out, "displacement.nii.gz"), True)
    sitk.WriteImage(sitk.CheckerBoard(sitk.RescaleIntensity(fixed, 0, 255),
                                      sitk.RescaleIntensity(moved, 0, 255), [8, 8, 1]),
                    os.path.join(a.out, "checker.nii.gz"), True)

    for p in a.also_warp:
        # label 必须用最近邻, 否则插值出小数标签
        is_label = "label" in os.path.basename(p).lower() or "mask" in os.path.basename(p).lower()
        img = sitk.ReadImage(p)
        r = sitk.Resample(img, fixed, tx,
                          sitk.sitkNearestNeighbor if is_label else sitk.sitkBSpline,
                          0.0, img.GetPixelID())
        out_p = os.path.join(a.out, "warped_" + os.path.basename(p))
        sitk.WriteImage(r, out_p, True)
        print(f"[ok] {out_p}")

    print(f"\n配准完成 -> {a.out}")


if __name__ == "__main__":
    run(build_parser(argparse.ArgumentParser()).parse_args())
