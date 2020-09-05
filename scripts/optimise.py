#!/usr/bin/env python

# Author: William Hunter
# This script handles both number of iterations and change stop criteria runs

# Import required modules:
from sys import argv
import topy
from topy import parser


def optimise(fname):
    # Set up ToPy:
    d = parser.tpd_file2dict(fname)
    if d['TO_TYPE'] == "trad":
        t = topy.TopologyTrad(topydict=d)
    elif d['TO_TYPE'] == "gen":
        t = topy.TopologyGen(topydict=d)
    else:
        raise ValueError("Unknown topology optimization approach: {}".format(d['TO_TYPE']))

    t.set_top_params()
    topy.optimise(t)


if __name__ == '__main__':
    optimise(argv[1])

