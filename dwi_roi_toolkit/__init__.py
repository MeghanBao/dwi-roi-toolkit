"""
dwi_roi_toolkit — semi-automatic whole-volume ROI delineation on diffusion-weighted MRI.

Modules
-------
convert     DICOM series -> NIfTI (b-value aware, PHI-stripping optional)
register    Rigid + B-spline deformable registration
candidates  Hyperintense focus detection (screening only, NOT diagnosis)
segment     Seed-based semi-automatic ROI growth -> ITK-SNAP label
"""

__version__ = "0.1.0"
__all__ = ["convert", "register", "candidates", "segment"]
