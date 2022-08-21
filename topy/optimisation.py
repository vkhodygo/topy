"""
# =============================================================================
# Optimization loop and topology preprocessing (for 'gen' method).
#
# Author: William Hunter, Tarcísio L. de Oliveira
# Copyright (C) 2008, 2015, William Hunter.
# Copyright (C) 2020, 2021, Tarcísio L. de Oliveira
# =============================================================================
"""

# -*- coding: utf-8 -*-
from os import path, makedirs
from time import time

from numpy import array
from mpi4py import MPI

from .utils import get_logger
from .visualisation import *

import cProfile
import pstats

logger = get_logger(__name__)


__all__ = ["optimise"]

def optimise(topology, save=True, dir='./iterations', apikey=''):
    # type: (Topology, bool, str) -> None
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        if not path.exists(dir):
            makedirs(dir)
        etas_avg = []

    # Optimising function:
    def _optimise_trad(t):
        t.fea()
    
        if rank == 0:
            t.sens_analysis()
            t.filter_sens_sigmund()
            t.update_desvars_oc()

            # Below this line we print info and create images or geometry:
            if t.nelz:
                params = {
                    'prefix': t.probname,
                    'iternum': t.itercount,
                    'time': 'none',
                    'dir': dir
                }
                if save:
                    create_3d_geom(t.desvars, **params)
            else:
                params = {
                    'prefix': t.probname,
                    'iternum': t.itercount,
                    'time': 'none',
                    'filetype': 'png',
                    'dir': dir
                }
                if save:
                    create_2d_imag(t.desvars, **params)

            
            if topology.topydict["PROB_TYPE"] != "heat":
                str_ = '%4i  | %3.6e | %3.3f | %3.4e | %3.3f | %3.3f |  %1.3f  |   %1.4f   | %3.3f'
                format_ = (t.itercount, t.objfval, t.desvars.mean(),\
                    t.change, t.p, t.q, t.eta.mean(), t.svtfrac, t.stress*1e-6)
                logger.info(str_ % format_)
            else:
                str_ = '%4i  | %3.6e | %3.3f | %3.4e | %3.3f | %3.3f |  %1.3f  |   %1.4f   '
                format_ = (t.itercount, t.objfval, t.desvars.mean(),\
                    t.change, t.p, t.q, t.eta.mean(), t.svtfrac)
                logger.info(str_ % format_)
            # Build a list of average etas:
            etas_avg.append(t.eta.mean())

    def _optimise_gen(t, Kfree):
        t.fea(Kfree)

        if rank == 0:
            t.sens_analysis()
            t.filter_sens_sigmund()
            t.update_desvars_oc()
            Kfree = t.updateK()

            # Below this line we print info and create images or geometry:
            if t.nelz:
                params = {
                    'prefix': t.probname,
                    'iternum': t.itercount,
                    'time': 'none',
                    'dir': dir
                }
                if save:
                    create_3d_geom(t.desvars, **params)
            else:
                params = {
                    'prefix': t.probname,
                    'iternum': t.itercount,
                    'time': 'none',
                    'filetype': 'png',
                    'dir': dir
                }
                if save:
                    create_2d_imag(t.desvars, **params)

            
            if topology.topydict["PROB_TYPE"] != "heat":
                str_ = '%4i  | %3.6e | %3.3f | %3.4e | %3.3f | %3.3f |  %1.3f  |   %1.4f   | %3.3f'
                format_ = (t.itercount, t.objfval, t.desvars.mean(),\
                    t.change, t.p, t.q, t.eta.mean(), t.svtfrac, t.stress*1e-6)
                logger.info(str_ % format_)
            else:
                str_ = '%4i  | %3.6e | %3.3f | %3.4e | %3.3f | %3.3f |  %1.3f  |   %1.4f   '
                format_ = (t.itercount, t.objfval, t.desvars.mean(),\
                    t.change, t.p, t.q, t.eta.mean(), t.svtfrac)
                logger.info(str_ % format_)
            # Build a list of average etas:
            etas_avg.append(t.eta.mean())

        return Kfree

    Kfree = None
    if rank == 0:
        ti = time()

    if topology.topydict["TO_TYPE"] == "gen":
        Kfree = topology.preprocess_space()
        if rank == 0:
            tp = time()

    if rank == 0:
        if topology.topydict["PROB_TYPE"] != "heat":
            # Create (plot) initial design domain:
            logger.info('\n' + '='*90)
            # Start optimisation runs, create rest of design domains:
            str_ = '%5s | %11s | %5s | %10s | %5s | %5s | %7s | %7s | %5s'
            format_ = ('Iter', 'Obj. func.  ', 'Vol. ', 'Change    ', \
                'P_FAC', 'Q_FAC', 'Ave ETA', 'Vol. frac.', 'Stress')
            logger.info(str_ % format_)
            logger.info('-'*90)
        else:
            # Create (plot) initial design domain:
            logger.info('\n' + '='*80)
            # Start optimisation runs, create rest of design domains:
            str_ = '%5s | %11s | %5s | %10s | %5s | %5s | %7s | %7s '
            format_ = ('Iter', 'Obj. func.  ', 'Vol. ', 'Change    ', \
                'P_FAC', 'Q_FAC', 'Ave ETA', 'Vol. frac.')
            logger.info(str_ % format_)
            logger.info('-'*80)

    # Optimize, and check for stop conditions
    if topology.topydict["TO_TYPE"] == "trad":
        if topology.chgstop > -1:
            while topology.change > topology.chgstop:
                _optimise_trad(topology)
        else:
            for i in range(topology.numiter):
                _optimise_trad(topology)
    else:
        if topology.chgstop > -1:
            while topology.change > topology.chgstop:
                Kfree = _optimise_gen(topology, Kfree)
        else:
            for i in range(topology.numiter):
                Kfree = _optimise_gen(topology, Kfree)

    if rank == 0:
        te = time()

        # Print iteration info:
        if topology.topydict['TO_TYPE'] == 'trad':
            logger.info('%d iterations took %3.3f minutes (%3.3f seconds/iteration)'\
                %(topology.itercount, (te - ti) / 60, (te - ti) / topology.itercount))
        else:
            logger.info('%d iterations took %3.3f minutes (%3.3f seconds/iteration, %3.3f seconds in preprocessing stage)'\
                %(topology.itercount, (te - ti) / 60, (te - tp) / topology.itercount, (tp - ti)))
        logger.info('Average of all ETA\'s = %3.3f (average of all a\'s = %3.3f)' \
            % (array(etas_avg).mean(), 1/array(etas_avg).mean() - 1))

