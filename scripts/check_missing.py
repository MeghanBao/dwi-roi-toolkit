#!/usr/bin/env python3
"""
收片对账 / Ingestion completeness check.

拿 DICOMDIR (刻盘/导出时生成的目录索引, 只有几 MB、不含像素, 但记录了本次检查
"应该"有哪些序列、每个序列多少张图) 跟实际文件夹里收到的图像逐序列对账, 一眼看出
缺哪些序列、哪些序列不全。

用法:
    python scripts/check_missing.py <DICOMDIR> <收到的文件夹>
    python scripts/check_missing.py <DICOMDIR> <收到的文件夹> --per-bvalue   # 逐 b 值展开

为什么用 DICOMDIR 而不是数文件:
    DICOMDIR 是发方导出时写的权威清单。文件夹里的图可能被截断、漏拷、或只导出了
    一层, 光数文件夹永远不知道"本该有多少"。对账才能发现层面被截断这类静默缺失。

坑 (已处理): 飞利浦会在每个序列里塞一个私有 series-data 文件 (SOP 名形如 `XX_*`,
    没有 Rows/Columns), 若把它当图像统计, 完全没到的序列会显示成"已收到 1 张", 看着
    像"不全"而不是"完全没到"。用 `hasattr(ds, "Rows")` 过滤掉这类非图像对象。

退出码: 有任何缺失/不全返回 1, 全部完整返回 0 (方便每收一例先跑一遍 / 接进 CI)。
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

import pydicom

# 与 convert.py 保持一致: 标准 b 值标签优先, 回退到飞利浦私有标签
STD_BVAL = (0x0018, 0x9087)
PHILIPS_BVAL = (0x2001, 0x1003)


def get_bval(ds):
    """取 b 值, 标准标签优先, 回退到飞利浦私有标签; 取不到返回 None"""
    for tag in (STD_BVAL, PHILIPS_BVAL):
        if tag in ds:
            try:
                v = ds[tag].value
                v = v[0] if isinstance(v, (list, tuple)) else v
                return int(round(float(v)))
            except Exception:
                pass
    return None


def read_dicomdir(dicomdir_path):
    """解析 DICOMDIR, 返回 {SeriesInstanceUID: {number, desc, modality, n_index}}"""
    ds = pydicom.dcmread(dicomdir_path, force=True)
    if "DirectoryRecordSequence" not in ds:
        raise SystemExit(f"{dicomdir_path} 不是有效的 DICOMDIR (没有 DirectoryRecordSequence)")

    index = {}
    current_uid = None
    for rec in ds.DirectoryRecordSequence:
        rtype = getattr(rec, "DirectoryRecordType", "")
        if rtype == "SERIES":
            uid = getattr(rec, "SeriesInstanceUID", None)
            if uid is None:
                current_uid = None
                continue
            current_uid = uid
            index.setdefault(uid, {
                "number": getattr(rec, "SeriesNumber", ""),
                # SERIES 记录里 SeriesDescription 是可选的; 缺就先留空, 之后用实收文件补
                "desc": getattr(rec, "SeriesDescription", "") or "",
                "modality": getattr(rec, "Modality", ""),
                "n_index": 0,
            })
        elif rtype == "IMAGE" and current_uid is not None:
            index[current_uid]["n_index"] += 1
    return index


def scan_received(data_dir):
    """
    递归读取实际收到的图像, 按 SeriesInstanceUID 分组。
    过滤掉飞利浦私有 series-data 等非图像对象 (无 Rows), 否则会把"没到"误报成"收到 1 张"。
    返回 {uid: {"number", "desc", "files": [(path, bval), ...]}}
    """
    got = defaultdict(lambda: {"number": "", "desc": "", "files": []})
    skipped_nonimage = 0
    for root, _, files in os.walk(data_dir):
        for fn in files:
            if fn.upper() == "DICOMDIR":
                continue
            path = os.path.join(root, fn)
            try:
                ds = pydicom.dcmread(path, stop_before_pixels=True, force=True)
            except Exception:
                continue
            if not hasattr(ds, "SeriesInstanceUID"):
                continue
            if not hasattr(ds, "Rows"):        # 非图像对象 (飞利浦 XX_* 私有序列数据)
                skipped_nonimage += 1
                continue
            uid = ds.SeriesInstanceUID
            entry = got[uid]
            entry["number"] = getattr(ds, "SeriesNumber", "")
            entry["desc"] = getattr(ds, "SeriesDescription", "") or entry["desc"]
            entry["files"].append((path, get_bval(ds)))
    return got, skipped_nonimage


def _w(s):
    """字符串的终端显示宽度 (CJK 按 2 计), 用于对齐"""
    return sum(2 if ord(c) > 0x2E7F else 1 for c in str(s))


def _pad(s, width):
    return str(s) + " " * max(0, width - _w(s))


def print_table(rows, headers):
    widths = [max(_w(h), *(_w(r[i]) for r in rows)) if rows else _w(h)
              for i, h in enumerate(headers)]
    line = "  ".join(_pad(h, w) for h, w in zip(headers, widths))
    print(line)
    for r in rows:
        print("  ".join(_pad(c, w) for c, w in zip(r, widths)))


def report(index, received, per_bvalue):
    rows = []
    total_index = total_got = 0
    incomplete = []

    # 按序列号排序; DICOMDIR 里有的序列全部列出
    for uid, meta in sorted(index.items(), key=lambda kv: str(kv[1]["number"])):
        n_index = meta["n_index"]
        got = received.get(uid, {"files": [], "desc": ""})
        n_got = len(got["files"])
        desc = meta["desc"] or got.get("desc", "")

        if n_got == 0:
            status = "缺失"
        elif n_got < n_index:
            status = f"不全(差{n_index - n_got})"
            incomplete.append(uid)
        elif n_got > n_index:
            status = f"多{n_got - n_index}"
        else:
            status = "完整"

        rows.append([meta["number"], desc, str(n_index), str(n_got), status])
        total_index += n_index
        total_got += n_got

    rows.append(["合计", "", str(total_index), str(total_got),
                 f"缺 {total_index - total_got} 张" if total_index > total_got else "完整"])

    print_table(rows, ["序列", "描述", "索引", "实收", "状态"])

    # 实收里出现、但 DICOMDIR 没索引的序列 (导出对不上, 也要提示)
    extra = [uid for uid in received if uid not in index]
    if extra:
        print("\n⚠ 以下序列在文件夹里存在, 但 DICOMDIR 未索引 (导出目录与实际不一致):")
        for uid in extra:
            g = received[uid]
            print(f"    序列 {g['number']} {g['desc']}  实收 {len(g['files'])} 张")

    # 逐 b 值展开 (默认对所有不全的 DWI 序列展开, 揪出层面截断)
    if per_bvalue and incomplete:
        for uid in incomplete:
            got = received.get(uid, {"files": []})
            bvals = [b for _, b in got["files"] if b is not None]
            if not bvals:
                continue
            meta = index[uid]
            counts = Counter(bvals)
            # 每个 b 值应有层数: 索引总数 / b 值种类数 (整除时可靠)
            n_b = len(counts)
            expect = meta["n_index"] // n_b if n_b and meta["n_index"] % n_b == 0 else None
            print(f"\n序列 {meta['number']} {meta['desc']} 逐 b 值 "
                  f"(索引 {meta['n_index']} = {n_b} 个 b 值 × "
                  f"{expect if expect else '?'} 层):")
            brows = []
            for b in sorted(counts):
                exp = str(expect) if expect else "?"
                brows.append([f"b{b}", str(counts[b]), exp,
                              "" if expect and counts[b] == expect else "← 截断"])
            print_table(brows, ["b 值", "实收", "应有", ""])

    return total_got < total_index or bool(extra)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="拿 DICOMDIR 索引对账实收文件, 找出缺失/不全的序列 "
                    "(Ingestion completeness check against DICOMDIR).")
    ap.add_argument("dicomdir", help="DICOMDIR 文件路径 (导出/刻盘时生成的目录索引)")
    ap.add_argument("data_dir", help="实际收到的 DICOM 文件夹")
    ap.add_argument("--per-bvalue", action="store_true",
                    help="对不全的 DWI 序列逐 b 值展开, 暴露层面截断")
    a = ap.parse_args(argv)

    index = read_dicomdir(a.dicomdir)
    received, skipped = scan_received(a.data_dir)
    if skipped:
        print(f"(已跳过 {skipped} 个非图像对象 / 飞利浦私有序列数据)\n")

    missing = report(index, received, a.per_bvalue)
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
