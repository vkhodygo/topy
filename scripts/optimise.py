#!/usr/bin/env python
## Author: William Hunter
"""
Optimze a TPD file from the commandline. See `python scripts/optimize.py --help`.

This script handles both number of iterations and change stop criteria runs.
"""

# Import required modules:
import topy
import argparse

parser = argparse.ArgumentParser(description="Optimize a TPD file.")
parser.add_argument("filename", type=str, help="the path to the TPD file")
parser.add_argument(
    "--vtk-format",
    default="binary",
    help='["binary"|"ascii"]: specify the format for VTK output files',
)


def optimise(filename, vtk_format="binary"):
    # type: (str, str) -> None
    """Optimize the TPD file at `fname`."""
    # Set up ToPy:
    t = topy.Topology()
    t.load_tpd_file(filename)
    t.set_top_params()
    topy.optimise(t)


if __name__ == "__main__":
    args = parser.parse_args()
    optimise(**args.__dict__)
