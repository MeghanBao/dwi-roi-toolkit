#!/usr/bin/env python3
"""
DICOM 序列 -> NIfTI (ITK-SNAP 可直接打开)

用法:
    python dcm2nii.py <dicom_dir> <out_dir> [--anonymize]

行为:
    - 递归扫描 dicom_dir, 按 SeriesInstanceUID 分组
    - DWI 序列按 b 值 (0018,9087 / Philips 私有 2001,1003) 再拆分成 b0/b500/b800...
    - 保留 spacing / origin / direction, 保证与 ITK-SNAP、后续配准坐标一致
    - 输出 <out_dir>/<SeriesNumber>_<描述>[_b<值>].nii.gz + 一份 series_index.csv
"""
import sys, os, re, csv, argparse
from collections import defaultdict
import numpy as np
import pydicom
import SimpleITK as sitk

# Philips 私有 b 值标签, 部分机器不写标准的 (0018,9087)
PHILIPS_BVAL = (0x2001, 0x1003)
STD_BVAL = (0x0018, 0x9087)

# 需要清除的 PHI 字段 (--anonymize 时生效)
PHI_TAGS = [
    "PatientName", "PatientID", "PatientBirthDate", "OtherPatientIDs",
    "PatientAddress", "PatientTelephoneNumbers", "InstitutionName",
    "InstitutionAddress", "ReferringPhysicianName", "PerformingPhysicianName",
    "OperatorsName", "AccessionNumber", "StudyID",
]


def safe(s):
    """把序列描述变成安全的文件名"""
    return re.sub(r"[^\w\-.]+", "_", str(s or "series")).strip("_")[:60]


def get_bval(ds):
    """取 b 值, 标准标签优先, 回退到 Philips 私有标签"""
    for tag in (STD_BVAL, PHILIPS_BVAL):
        if tag in ds:
            try:
                v = ds[tag].value
                v = v[0] if isinstance(v, (list, tuple)) else v
                return int(round(float(v)))
            except Exception:
                pass
    return None


def scan(dicom_dir):
    """递归读取所有 DICOM, 按 (SeriesInstanceUID, b值) 分组"""
    groups = defaultdict(list)
    for root, _, files in os.walk(dicom_dir):
        for fn in files:
            path = os.path.join(root, fn)
            try:
                # stop_before_pixels: 只读头, 分组阶段不需要像素, 快很多
                ds = pydicom.dcmread(path, stop_before_pixels=True, force=True)
            except Exception:
                continue
            if not hasattr(ds, "SeriesInstanceUID"):
                continue
            groups[(ds.SeriesInstanceUID, get_bval(ds))].append(path)
    return groups


def sort_slices(paths):
    """按 ImagePositionPatient 沿层面法向排序; 缺失时回退 InstanceNumber"""
    metas = []
    for p in paths:
        ds = pydicom.dcmread(p, stop_before_pixels=True, force=True)
        metas.append((p, ds))
    try:
        iop = np.array(metas[0][1].ImageOrientationPatient, dtype=float)
        normal = np.cross(iop[:3], iop[3:])  # 层面法向量
        metas.sort(key=lambda m: float(np.dot(np.array(m[1].ImagePositionPatient, float), normal)))
    except Exception:
        metas.sort(key=lambda m: int(getattr(m[1], "InstanceNumber", 0)))
    return [p for p, _ in metas]


def anonymize_inplace(reader_files, out_dir):
    """写一份去标识的 DICOM 副本(可选), 原始文件不改动"""
    anon_dir = os.path.join(out_dir, "dicom_anon")
    os.makedirs(anon_dir, exist_ok=True)
    for i, p in enumerate(reader_files):
        ds = pydicom.dcmread(p, force=True)
        for t in PHI_TAGS:
            if hasattr(ds, t):
                setattr(ds, t, "")
        ds.PatientIdentityRemoved = "YES"
        ds.save_as(os.path.join(anon_dir, f"{i:05d}.dcm"))
    return anon_dir


def convert(dicom_dir, out_dir, anonymize=False):
    os.makedirs(out_dir, exist_ok=True)
    groups = scan(dicom_dir)
    if not groups:
        raise SystemExit(f"在 {dicom_dir} 下没找到 DICOM 文件")

    index_rows = []
    for (uid, bval), paths in sorted(groups.items(), key=lambda kv: (kv[0][0], -1 if kv[0][1] is None else kv[0][1])):
        paths = sort_slices(paths)
        ref = pydicom.dcmread(paths[0], stop_before_pixels=True, force=True)

        # SimpleITK 的 ImageSeriesReader 会正确处理 spacing / direction / rescale
        reader = sitk.ImageSeriesReader()
        reader.SetFileNames(paths)
        try:
            img = reader.Execute()
        except RuntimeError as e:
            # 非图像 DICOM(二次截图/演示状态等)无法读取; 跳过而非整体中止
            print(f"[skip] {getattr(ref,'SeriesDescription','')!r} "
                  f"(uid={uid[-12:]}, b={bval}): 无法作为图像读取 ({str(e).splitlines()[-1]})")
            continue

        # 单层序列 SimpleITK 会给出 z-spacing=1, 用 DICOM 的层间距修正
        if img.GetSize()[2] == 1:
            sp = list(img.GetSpacing())
            sp[2] = float(getattr(ref, "SpacingBetweenSlices", None)
                          or getattr(ref, "SliceThickness", 1) or 1)
            img.SetSpacing(sp)

        name = f"{getattr(ref,'SeriesNumber','NA')}_{safe(getattr(ref,'SeriesDescription',''))}"
        if bval is not None:
            name += f"_b{bval}"
        out_path = os.path.join(out_dir, name + ".nii.gz")
        sitk.WriteImage(img, out_path, useCompression=True)

        index_rows.append({
            "file": os.path.basename(out_path),
            "series_number": getattr(ref, "SeriesNumber", ""),
            "description": getattr(ref, "SeriesDescription", ""),
            "b_value": bval if bval is not None else "",
            "n_slices": img.GetSize()[2],
            "size": "x".join(map(str, img.GetSize())),
            "spacing": "x".join(f"{s:.4f}" for s in img.GetSpacing()),
            "TR": getattr(ref, "RepetitionTime", ""),
            "TE": getattr(ref, "EchoTime", ""),
        })
        print(f"[ok] {out_path}  size={img.GetSize()}  spacing={tuple(round(s,3) for s in img.GetSpacing())}")

        if anonymize:
            anonymize_inplace(paths, out_dir)

    with open(os.path.join(out_dir, "series_index.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(index_rows[0].keys()))
        w.writeheader()
        w.writerows(index_rows)
    print(f"\n共 {len(index_rows)} 个序列 -> {out_dir}/series_index.csv")


def build_parser(ap):
    ap.add_argument("dicom_dir", help="directory containing the DICOM series")
    ap.add_argument("out_dir", help="output directory for NIfTI files")
    ap.add_argument("--anonymize", action="store_true",
                    help="also write a PHI-stripped copy of the DICOM files")
    return ap


def run(a):
    convert(a.dicom_dir, a.out_dir, a.anonymize)


if __name__ == "__main__":
    run(build_parser(argparse.ArgumentParser()).parse_args())
