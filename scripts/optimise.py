#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: William Hunter
# This script handles both number of iterations and change stop criteria runs

# Import required modules:
import signal
from sys import argv
import topy



def optimise(fname, outdir, apikey):

    # Set up ToPy:
    t = topy.Topology()
    t.load_tpd_file(fname)
    t.set_top_params()
    topy.optimise(t, dir=outdir, apikey=apikey)


if __name__ == '__main__':

    optimise(argv[1], argv[2], argv[3])

