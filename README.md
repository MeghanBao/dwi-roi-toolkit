# dwi-roi-toolkit

[English](README.md) | [简体中文](README.zh-CN.md)

Semi-automatic whole-volume ROI delineation and quantification on diffusion-weighted MRI.

Turns lesion tracing from a slice-by-slice manual chore into two steps: a radiologist places a
seed point, the toolkit grows it into a 3D ROI, and the result opens directly in ITK-SNAP for
review and correction.

> **Research use only.** Not a diagnostic device, no regulatory clearance. The toolkit does not
> decide what is or is not a lesion — a qualified radiologist supplies the target and must review
> every output before use.

---

## Why

Prospective DWI studies need per-case volume and ADC measurements. Producing them by hand means
tracing the lesion on every slice it spans, for every case in the cohort. The bottleneck is not
judgment — a radiologist locates the lesion in seconds — it is the tracing.

This toolkit keeps the judgment with the clinician and automates the tracing:

```
DICOM in  →  deformable registration  →  seed point (clinician)  →  3D ROI  →  volume + ADC out
```

---

## Install

```bash
git clone https://github.com/<you>/dwi-roi-toolkit.git
cd dwi-roi-toolkit
pip install -e .
```

Python ≥ 3.9. Dependencies: SimpleITK, pydicom, numpy, scipy.

Installs a `dwi-roi` command; `python -m dwi_roi_toolkit` works equivalently.

---

## Pipeline

### 1. Convert

```bash
dwi-roi convert ./patient_dicom ./case01/nifti --anonymize
```

Recursively scans for DICOM files, groups by `SeriesInstanceUID`, then splits diffusion series by
b-value — reading the standard tag `(0018,9087)` and falling back to the Philips private tag
`(2001,1003)`, which some scanners populate instead. Slices are ordered by projecting
`ImagePositionPatient` onto the slice normal rather than trusting `InstanceNumber`.

Spacing, origin and direction are preserved, so every downstream step and ITK-SNAP itself share
one physical coordinate system.

`--anonymize` also writes a PHI-stripped copy of the DICOM files.

### 2. Register

```bash
dwi-roi register --fixed T2.nii.gz --moving DWI_b1000.nii.gz --out ./case01/reg \
                 --also-warp ADC.nii.gz
```

Thoracic DWI and structural images are acquired at different points in the respiratory and cardiac
cycle, so rigid alignment alone leaves lung parenchyma several millimetres off. The pipeline runs
centroid initialisation → `VersorRigid3D` → B-spline free-form deformation, each over a three-level
resolution pyramid with Mattes mutual information (DWI vs. T2 is multi-modal; sum-of-squares would
not converge).

Control points sit roughly every 30 mm — coarse enough not to fit noise, fine enough to absorb
respiratory motion.

Outputs `moved.nii.gz`, a reusable `transform.tfm`, the displacement field, and a checkerboard
image for visual QC. `--also-warp` applies the same transform to ADC maps, other b-values or
existing labels (labels get nearest-neighbour resampling automatically).

`--mode demons` offers diffeomorphic Demons for single-modality cases.

### 3. Detect candidates

```bash
dwi-roi candidates --image ./case01/reg/moved.nii.gz --out ./case01 --topk 8
```

Builds a body mask, smooths to suppress EPI noise spikes, thresholds at a high intra-body
percentile, and reports connected components with peak-intensity coordinates.

**This is a screening aid, not a detector.** Normal structures — spinal cord, vessels, bowel — are
routinely the brightest foci on high-b-value DWI. The output is a shortlist for a clinician to
choose from, written to `candidates.csv`.

### 4. Segment

```bash
dwi-roi segment --image ./case01/reg/moved.nii.gz --seed 223 212 44 \
                --adc ./case01/reg/warped_ADC.nii.gz --out ./case01/roi
```

Curvature anisotropic diffusion denoising → `ConfidenceConnected` region growing → threshold
level-set boundary refinement → morphological cleanup → retain the seed's connected component.

Region growing adapts its threshold from seed-neighbourhood statistics rather than using a fixed
intensity cutoff, which matters because DWI intensity is not standardised across scanners or
sessions.

Writes `label.nii.gz` (uint8, values in {0,1}) and `roi_stats.json` with volume in mm³ and mL,
signal statistics, and — when `--adc` is given — ADCmean, ADCsd and ADCmin over the ROI.

If the level set collapses the contour or lets it explode, the result falls back to the region-
growing mask with a warning rather than silently emitting an empty or blown-out ROI.

**Tuning:** `--multiplier` controls how permissively the region grows; lower it if the ROI leaks
into surrounding tissue. `--propagation` and `--curvature` control the level-set forces.

### 5. Review in ITK-SNAP

```
File > Open Main Image            →  image.nii.gz
Segmentation > Open Segmentation  →  label.nii.gz
```

The label shares size, spacing, origin and direction with the image, so it overlays exactly. Brush,
eraser and 3D ball tools all work; `Segmentation > Save Segmentation` writes the corrected version
back.

---

## Batch processing

Everything except seed selection is scriptable:

```bash
for case in cohort/*/; do
  dwi-roi convert  "$case/dicom" "$case/nifti"
  dwi-roi register --fixed "$case/nifti/T2.nii.gz" \
                   --moving "$case/nifti/DWI_b1000.nii.gz" \
                   --out "$case/reg" --also-warp "$case/nifti/ADC.nii.gz"
  dwi-roi candidates --image "$case/reg/moved.nii.gz" --out "$case"
done
# clinician records seed coordinates, then:
dwi-roi segment --image "$case/reg/moved.nii.gz" --seed $X $Y $Z \
                --adc "$case/reg/warped_ADC.nii.gz" --out "$case/roi"
```

---

## Input requirements

| Input | Needed for |
|---|---|
| Full DWI series, all b-values | Whole-volume ROI; a single slice yields only a 2D contour |
| ADC map (or ≥ 2 b-values) | ADC statistics — usually the primary quantitative endpoint |
| Structural series (T2 or CT) | Deformable registration target |

A common failure mode is receiving one exported slice instead of the study directory. Whole-volume
delineation is impossible in that case.

### Reconcile against DICOMDIR before ingesting

An exported study ships with a `DICOMDIR` — a few-megabyte index (no pixels) listing every series
the study *should* contain and how many images each holds. Counting files in the folder never tells
you what was *supposed* to arrive, so silent truncation slips through. Reconcile the two:

```bash
python scripts/check_missing.py /path/to/DICOMDIR /path/to/received_folder --per-bvalue
```

It prints an indexed-vs-received table per series (`缺失` / `不全(差N)` / `完整`) and, with
`--per-bvalue`, expands incomplete diffusion series by b-value — which is how a DKI series that
looked like it had merely "fewer slices" turned out to be truncated (each b-value short 7–8 slices).
Non-image objects such as Philips private series-data files (`XX_*`, no `Rows`) are filtered out, so
a series that never arrived reads as `缺失`, not a misleading "received 1". The script exits non-zero
when anything is missing, so it drops straight into an intake step. Run it on every new case.

---

## Patient data and this repository

This repository contains **no patient data** — no DICOM files, no NIfTI volumes, no rendered images
of any scan.

Medical images remain protected health information after DICOM tags are stripped, and git history
cannot be reliably scrubbed once pushed. `.gitignore` blocks medical image formats, common data
directory names, image renders and archives. A stricter guard is available as a pre-commit hook:

```bash
ln -s ../../scripts/check_no_phi.py .git/hooks/pre-commit
chmod +x scripts/check_no_phi.py
python scripts/check_no_phi.py --all   # scan the whole tree
```

It rejects blocked extensions and sniffs for the `DICM` magic header, since DICOM files are
frequently distributed without any extension.

Use `dwi-roi convert --anonymize` before sharing data with collaborators.

---

## Limitations

These are load-bearing, not boilerplate.

- **Segmentation is parameter-sensitive.** On a noisy single-slice test, moving the seed by one
  voxel changed the resulting volume from empty to several thousand mm³. Full volumes are better
  conditioned — region growing then estimates its statistics from a genuine 3D neighbourhood — but
  every cohort needs its parameters fixed up front and held constant.
- **Detection ranks bright voxels and nothing more.** It has no notion of anatomy and will happily
  rank the spinal cord first. It also cannot tell you that a slice contains no lesion at all; asked
  for candidates on a normal slice, it returns normal structures.
- **Registration is tuned for thoracic anatomy.** Other body regions need different control-point
  spacing.
- **No validation against expert manual segmentation.** Dice agreement has not been measured. Doing
  so is a prerequisite for any study relying on these outputs, and no accuracy claim should be made
  until it is.

---

## License

MIT — see [LICENSE](LICENSE).
