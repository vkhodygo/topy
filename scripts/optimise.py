#!/usr/bin/env python

# Author: William Hunter, Tarcísio L. de Oliveira
# This script handles both number of iterations and change stop criteria runs

# Import required modules:
from sys import argv
import topy
from topy import parser
from mpi4py import MPI
from mumps import DMumpsContext


def optimise(fname):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # MPI hack to prevent lambda serialization (unsupported by mpi4py)
    t = topy.TopologyTrad()
    t.topydict['TO_TYPE'] = "trad"
    t.chgstop = -1
    t.numiter = 0
    t.probtype = "comp"
    gen = True
    # Set up ToPy:
    if rank == 0:
        d = parser.tpd_file2dict(fname)

        if d['TO_TYPE'] == "trad":
            t = topy.TopologyTrad(topydict=d)
            gen = False
        elif d['TO_TYPE'] == "gen":
            t = topy.TopologyGen(topydict=d)
        else:
            raise ValueError("Unknown topology optimization approach: {}".format(d['TO_TYPE']))

        t.set_top_params()
        t.chgstop = -1

    gen = comm.bcast(gen, root=0)
    if gen and rank != 0:
        t = topy.TopologyGen()
        t.topydict['TO_TYPE'] = "gen"
        t.probtype = "comp"
        t.chgstop = -1
        t.numiter = 0
    t.topydict['TO_TYPE'] = comm.bcast(t.topydict['TO_TYPE'], root=0)
    t.probtype = comm.bcast(t.probtype, root=0)
    t.chgstop = comm.bcast(t.chgstop, root=0)
    t.numiter = comm.bcast(t.numiter, root=0)

    topy.optimise(t)
    
if __name__ == '__main__':
    optimise(argv[1])

