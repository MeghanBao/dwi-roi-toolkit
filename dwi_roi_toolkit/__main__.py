"""Unified CLI: python -m dwi_roi_toolkit <subcommand> ..."""
import argparse
import sys

from . import __version__, candidates, convert, register, segment

SUBCOMMANDS = {
    "convert": (convert, "DICOM series -> NIfTI (splits by series and b-value)"),
    "register": (register, "Rigid + B-spline deformable registration"),
    "candidates": (candidates, "Detect hyperintense foci (screening aid, not diagnosis)"),
    "segment": (segment, "Seed-based semi-automatic ROI -> ITK-SNAP label"),
}


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="dwi-roi",
        description="Semi-automatic whole-volume ROI delineation and quantification on diffusion-weighted MRI. "
                    "Research use only - not a diagnostic device.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Typical pipeline:\n"
               "  dwi-roi convert    ./dicom ./case01/nifti\n"
               "  dwi-roi register   --fixed T2.nii.gz --moving DWI_b1000.nii.gz --out ./case01/reg\n"
               "  dwi-roi candidates --image ./case01/reg/moved.nii.gz --out ./case01\n"
               "  # radiologist picks the target seed from candidates.csv\n"
               "  dwi-roi segment    --image ./case01/reg/moved.nii.gz --seed X Y Z --out ./case01/roi\n",
    )
    ap.add_argument("--version", action="version", version=f"dwi-roi-toolkit {__version__}")
    sub = ap.add_subparsers(dest="cmd", metavar="<subcommand>")

    for name, (mod, help_text) in SUBCOMMANDS.items():
        p = sub.add_parser(name, help=help_text, description=help_text)
        mod.build_parser(p)
        p.set_defaults(_run=mod.run)

    args = ap.parse_args(argv)
    if not getattr(args, "cmd", None):
        ap.print_help()
        return 1
    args._run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
