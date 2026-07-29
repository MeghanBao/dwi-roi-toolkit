#!/usr/bin/env python3
"""
DWI 高信号候选灶检出 —— 供医生挑种子点用, 不做诊断判定。

用法:
    python find_candidates.py --image DWI_b1000.nii.gz --out out_dir [--topk 8]

原理: 体部掩膜内做 LoG 斑点增强 + 自适应阈值, 列出高信号连通域及其体素坐标,
医生在标注图上认哪个是靶病灶, 把对应坐标喂给 segment_roi.py --seed。
"""
import os, csv, argparse
import numpy as np
import SimpleITK as sitk
from scipy import ndimage


def body_mask(arr):
    """体部掩膜: 排除背景空气与体外噪声"""
    thr = np.percentile(arr, 40)
    m = arr > thr
    m = ndimage.binary_closing(m, np.ones((1, 5, 5)))
    m = ndimage.binary_fill_holes(m)
    lab, n = ndimage.label(m)
    if n > 1:                                    # 只留最大连通域(躯干)
        sizes = ndimage.sum(m, lab, range(1, n + 1))
        m = lab == (np.argmax(sizes) + 1)
    return m


def build_parser(ap):
    ap.add_argument("--image", required=True, help="DWI image (NIfTI), preferably high b-value")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--topk", type=int, default=8, help="max number of candidates to report")
    ap.add_argument("--percentile", type=float, default=99.0,
                    help="intensity percentile (within body mask) used as threshold")
    ap.add_argument("--min-mm3", type=float, default=20.0, help="minimum candidate volume")
    return ap


def run(a):
    os.makedirs(a.out, exist_ok=True)

    img = sitk.ReadImage(a.image, sitk.sitkFloat32)
    arr = sitk.GetArrayFromImage(img)            # [z, y, x]
    sp = img.GetSpacing()                        # (x, y, z)
    vox = float(np.prod(sp))

    body = body_mask(arr)
    # 轻度高斯平滑抑制 EPI 噪声尖峰, 再取体内高分位阈值
    sm = ndimage.gaussian_filter(arr, sigma=(0, 1.0, 1.0))
    thr = np.percentile(sm[body], a.percentile)
    hi = (sm > thr) & body

    lab, n = ndimage.label(hi)
    rows = []
    for i in range(1, n + 1):
        m = lab == i
        v = m.sum() * vox
        if v < a.min_mm3:
            continue
        zs, ys, xs = np.where(m)
        # 峰值点比质心更适合当区域生长的种子
        idx = np.argmax(arr[m])
        rows.append({
            "id": len(rows) + 1,
            "seed_x": int(xs[idx]), "seed_y": int(ys[idx]), "seed_z": int(zs[idx]),
            "centroid_x": int(round(xs.mean())), "centroid_y": int(round(ys.mean())),
            "centroid_z": int(round(zs.mean())),
            "voxels": int(m.sum()),
            "volume_mm3": round(v, 2),
            "max_signal": round(float(arr[m].max()), 1),
            "mean_signal": round(float(arr[m].mean()), 1),
        })
    rows.sort(key=lambda r: -r["max_signal"])
    rows = rows[:a.topk]
    for k, r in enumerate(rows, 1):
        r["id"] = k

    csv_path = os.path.join(a.out, "candidates.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["id"])
        w.writeheader(); w.writerows(rows)
    print(f"检出 {len(rows)} 个候选高信号区 -> {csv_path}")
    for r in rows:
        print(f"  #{r['id']}  seed=({r['seed_x']},{r['seed_y']},{r['seed_z']})  "
              f"{r['volume_mm3']}mm3  peak={r['max_signal']}")
    return rows


if __name__ == "__main__":
    run(build_parser(argparse.ArgumentParser()).parse_args())
