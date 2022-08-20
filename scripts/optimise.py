#!/usr/bin/env python
# -*- coding: utf-8 -*-

This script handles both number of iterations and change stop criteria runs.
"""

# Import required modules:
import signal
from sys import argv
import topy
import argparse

parser = argparse.ArgumentParser(description="Optimize a TPD file.")
parser.add_argument("filename", type=str, help="the path to the TPD file")
parser.add_argument(
    "--vtk-format",
    default="binary",
    help='["binary"|"ascii"]: specify the format for VTK output files',
)


def optimise(fname, outdir, apikey):

    # Set up ToPy:
    t = topy.Topology()
    t.load_tpd_file(filename)
    t.set_top_params()
    topy.optimise(t, dir=outdir, apikey=apikey)


if __name__ == '__main__':

    optimise(argv[1], argv[2], argv[3])

