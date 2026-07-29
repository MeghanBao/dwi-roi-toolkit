# dwi-roi-toolkit

[English](README.md) | [简体中文](README.zh-CN.md)

弥散加权 MRI 的半自动全体积 ROI 勾画与定量工具。

把病灶勾画从逐层手动的苦差事简化为两步：由放射科医生放置一个种子点，工具箱将其生长为
三维 ROI，结果直接在 ITK-SNAP 中打开以供审阅和修正。

> **仅供科研使用。** 非诊断设备，未获任何监管许可。本工具箱不判断什么是或不是病灶——由具备
> 资质的放射科医生提供目标，并且必须在使用前审阅每一个输出结果。

---

## 为什么

前瞻性 DWI 研究需要逐例的体积和 ADC 测量。手工完成意味着要在病灶跨越的每一层上进行勾画，
且队列中的每一例都要如此。瓶颈不在于判断——放射科医生几秒钟就能定位病灶——而在于勾画本身。

本工具箱把判断留给临床医生，并将勾画自动化：

```
DICOM 输入  →  形变配准  →  种子点（临床医生）  →  三维 ROI  →  输出体积 + ADC
```

---

## 安装

```bash
git clone https://github.com/<you>/dwi-roi-toolkit.git
cd dwi-roi-toolkit
pip install -e .
```

Python ≥ 3.9。依赖：SimpleITK、pydicom、numpy、scipy。

安装后提供 `dwi-roi` 命令；`python -m dwi_roi_toolkit` 效果等同。

---

## 处理流程

### 1. 转换

```bash
dwi-roi convert ./patient_dicom ./case01/nifti --anonymize
```

递归扫描 DICOM 文件，按 `SeriesInstanceUID` 分组，再按 b 值拆分弥散序列——读取标准标签
`(0018,9087)`，并回退到部分扫描仪改用的飞利浦私有标签 `(2001,1003)`。层面顺序按
`ImagePositionPatient` 在层法向上的投影排序，而非依赖 `InstanceNumber`。

间距、原点和方向都被保留，因此每一个下游步骤乃至 ITK-SNAP 本身，都共享同一个物理坐标系。

`--anonymize` 还会额外写出一份去除 PHI（受保护健康信息）的 DICOM 文件副本。

### 2. 配准

```bash
dwi-roi register --fixed T2.nii.gz --moving DWI_b1000.nii.gz --out ./case01/reg \
                 --also-warp ADC.nii.gz
```

胸部 DWI 与结构像是在呼吸与心动周期的不同时相采集的，因此仅靠刚性对齐会使肺实质偏差数毫米。
流程依次运行：质心初始化 → `VersorRigid3D` → B 样条自由形变，每一步都在三层分辨率金字塔上、
使用 Mattes 互信息进行（DWI 与 T2 属多模态；平方和度量无法收敛）。

控制点约每 30 mm 一个——足够粗以避免拟合噪声，又足够细以吸收呼吸运动。

输出 `moved.nii.gz`、一个可复用的 `transform.tfm`、位移场，以及用于目视质控的棋盘格图像。
`--also-warp` 会将同一变换应用于 ADC 图、其他 b 值或已有的标签（标签会自动采用最近邻重采样）。

`--mode demons` 为单模态场景提供微分同胚 Demons 配准。

### 3. 检测候选

```bash
dwi-roi candidates --image ./case01/reg/moved.nii.gz --out ./case01 --topk 8
```

构建体部掩膜，平滑以抑制 EPI 噪声尖峰，在体内高百分位处做阈值分割，并报告各连通分量及其
峰值强度坐标。

**这是筛查辅助，而非检测器。** 正常结构——脊髓、血管、肠道——在高 b 值 DWI 上常常是最亮的
病灶点。其输出是一份供临床医生挑选的候选清单，写入 `candidates.csv`。

### 4. 分割

```bash
dwi-roi segment --image ./case01/reg/moved.nii.gz --seed 223 212 44 \
                --adc ./case01/reg/warped_ADC.nii.gz --out ./case01/roi
```

曲率各向异性扩散去噪 → `ConfidenceConnected` 区域生长 → 阈值水平集边界细化 → 形态学清理 →
保留种子所在的连通分量。

区域生长从种子邻域统计量自适应地确定阈值，而非采用固定强度截断——这一点很关键，因为 DWI
强度在不同扫描仪或不同扫描间并不标准化。

写出 `label.nii.gz`（uint8，取值 {0,1}）和 `roi_stats.json`，其中包含以 mm³ 和 mL 表示的体积、
信号统计量，以及——当给出 `--adc` 时——ROI 内的 ADCmean、ADCsd 和 ADCmin。

如果水平集使轮廓塌缩或爆炸式扩张，结果会带警告回退到区域生长掩膜，而不是悄悄地输出一个空的
或膨胀失控的 ROI。

**调参：** `--multiplier` 控制区域生长的宽松程度；若 ROI 渗漏到周围组织，请调低它。
`--propagation` 和 `--curvature` 控制水平集的作用力。

### 5. 在 ITK-SNAP 中审阅

```
File > Open Main Image            →  image.nii.gz
Segmentation > Open Segmentation  →  label.nii.gz
```

标签与图像共享尺寸、间距、原点和方向，因此能精确叠加。画笔、橡皮和三维球体工具均可使用；
`Segmentation > Save Segmentation` 会把修正后的版本写回。

---

## 批处理

除种子选择外的一切都可脚本化：

```bash
for case in cohort/*/; do
  dwi-roi convert  "$case/dicom" "$case/nifti"
  dwi-roi register --fixed "$case/nifti/T2.nii.gz" \
                   --moving "$case/nifti/DWI_b1000.nii.gz" \
                   --out "$case/reg" --also-warp "$case/nifti/ADC.nii.gz"
  dwi-roi candidates --image "$case/reg/moved.nii.gz" --out "$case"
done
# 临床医生记录种子坐标，然后：
dwi-roi segment --image "$case/reg/moved.nii.gz" --seed $X $Y $Z \
                --adc "$case/reg/warped_ADC.nii.gz" --out "$case/roi"
```

---

## 输入要求

| 输入 | 用途 |
|---|---|
| 完整 DWI 序列，全部 b 值 | 全体积 ROI；单层只能得到二维轮廓 |
| ADC 图（或 ≥ 2 个 b 值） | ADC 统计量——通常是主要定量终点 |
| 结构序列（T2 或 CT） | 形变配准的目标图像 |

一种常见的失败情形是收到一张导出的单层图像，而非整个检查目录。此时无法进行全体积勾画。

---

## 患者数据与本仓库

本仓库**不含任何患者数据**——没有 DICOM 文件、没有 NIfTI 体数据、也没有任何扫描的渲染图像。

医学影像在剥离 DICOM 标签后仍属于受保护健康信息，而 git 历史一旦推送便无法可靠地清除。
`.gitignore` 拦截医学影像格式、常见数据目录名、图像渲染文件和压缩包。更严格的防护以
pre-commit 钩子形式提供：

```bash
ln -s ../../scripts/check_no_phi.py .git/hooks/pre-commit
chmod +x scripts/check_no_phi.py
python scripts/check_no_phi.py --all   # 扫描整个目录树
```

它会拒绝被封禁的扩展名，并嗅探 `DICM` 魔术头，因为 DICOM 文件常常在没有任何扩展名的情况下
分发。

在与协作者共享数据前，请使用 `dwi-roi convert --anonymize`。

---

## 局限性

以下都是实质性的、而非套话。

- **分割对参数敏感。** 在一个有噪声的单层测试中，把种子移动一个体素，就会使所得体积从空
  变为数千 mm³。全体积的条件更好——区域生长此时能从真实的三维邻域估计统计量——但每个队列都
  需要事先固定参数并保持恒定。
- **检测只是对亮体素排序，仅此而已。** 它没有解剖学概念，会毫不犹豫地把脊髓排在第一位。它也
  无法告诉你某一层根本不含病灶；若在正常层面上请求候选，它会返回正常结构。
- **配准针对胸部解剖调优。** 其他身体区域需要不同的控制点间距。
- **未与专家手工分割做验证。** Dice 一致性尚未测量。在任何依赖这些输出的研究中，这样做都是
  前提；在完成之前，不应做出任何准确性声明。

---

## 许可证

MIT——见 [LICENSE](LICENSE)。
