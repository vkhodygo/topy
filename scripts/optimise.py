#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: William Hunter
# This script handles both number of iterations and change stop criteria runs

# Import required modules:
import signal
from sys import argv
import topy


def handler(signum, frame):
     print("Forever is over!")
     #raise Exception("end of time")

def optimise(fname):
    # Set up ToPy:
    t = topy.Topology()
    t.load_tpd_file(fname)
    t.set_top_params()
    topy.optimise(t)


if __name__ == '__main__':
    profile = True
    if profile:
        #signal.signal(signal.SIGALRM, handler)
        #signal.alarm(100)

        import cProfile, pstats, io
        pr = cProfile.Profile()
        pr.enable()

        #try:
        optimise(argv[1])
        #except Exception:
            #pass

        pr.disable()
        s = io.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
    else:
        optimise(argv[1])

