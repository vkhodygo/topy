"""
# =============================================================================
# A class to optimise the topology of a design domain for defined boundary
# conditions. Data is read from an input file, see 'examples' directory.
# Uses traditional SIMP topology optimization.
#
# Authors: William Hunter, Tarcísio L. de Oliveira
# Copyright (C) 2008, 2015, William Hunter.
# Copyright (C) 2020, 2021, Tarcísio L. de Oliveira
# =============================================================================
"""
import os

import numpy as np
import numexpr as ne

from scipy import sparse
from scipy.sparse import linalg
from sympy import symbols
from mumps import DMumpsContext
from mpi4py import MPI

from .utils import get_logger
from .parser import tpd_file2dict, config2dict
from .visualisation import *

logger = get_logger(__name__)
__all__ = ['TopologyTrad']


MAX_ITERS = 250

SOLID, VOID = 1.000, 0.001 #  Upper and lower bound value for design variables
KDATUM = 0.1 #  Reference stiffness value of springs for mechanism synthesis

# Constants for exponential approximation:
A_LOW = -3 #  Lower restriction on 'a' for exponential approximation
A_UPP = -1e-5 #  Upper restriction on 'a' for exponential approximation


# =======================
# === ToPy base class ===
# =======================
class TopologyTrad:
    """
    A class to optimise the topology of a design domain for defined boundary
    values. Data is read from an input file (see 'examples' folder).

    """
    def __init__(self, config=None, topydict={}, pcount=0, 
                 qcount=0, itercount=0, change=1, svtfrac=None):
        self.pcount = pcount #  Counter for continuation of p
        self.qcount = qcount #  Counter for continuation of q for GSF
        self.itercount = itercount #  Internal counter
        self.change = change
        self.svtfrac = svtfrac

        if config:
            self.topydict = config2dict(config.copy())
        else:
            self.topydict = topydict #  Store tpd file data in dictionary


    # ======================
    # === Public methods ===
    # ======================
    def load_tpd_file(self, fname):
        """
        DEPRECATED. Need to check TO_TYPE before instantiating.
        Load a ToPy problem definition (TPD) file, return a dictionary:

        INPUTS:
            fname -- file name of tpd file.

        EXAMPLES:
            >>> import topy
            >>> t = topy.Topology()
            >>> # load a ToPy problem definition file:
            >>> t.load_tpd_file('filename.tpd')

        """
        self.tpdfname = fname
        self.topydict = tpd_file2dict(fname)

    def set_top_params(self):
        """
        Set topology optimisation problem parameters (you must already have
        instantiated a Topology object).

        EXAMPLES:
            >>> import topy
            >>> t = topy.Topology()
            >>> # load a ToPy problem definition file:
            >>> t.load_tpd_file('filename.tpd')
            >>> # now set the problem parameters:
            >>> t.set_top_params()

        You can now access all the class methods or instance variables. To see
        the dictionary of parameters, just type

            >>> t.topydict

        You can also reset some parameters without reloading the file, but use
        with caution. For example, you can change the filer radius (FILT_RAD)
        as follows:

            >>> t.topydict['FILT_RAD'] = 1.5

        Remember to reset the problem parameters again, whereafter you can
        solve the problem.

            >>> t.set_top_params()

        See also: load_tpd_file

        """
        # Set all the mandatory minimum amount of parameters that constitutes
        # a completely defined topology optimisation problem:
        if not self.topydict:
            raise Exception('You must first load a TPD file!')
        self.probtype = self.topydict['PROB_TYPE'] #  Problem type
        self.probname = self.topydict.get('PROB_NAME', '') #  Problem name
        self.volfrac = self.topydict['VOL_FRAC'] #  Volume fraction
        self.filtrad = self.topydict['FILT_RAD'] #  Filter radius
        self.p = self.topydict['P_FAC'] #  'Standard' penalisation factor
        self.dofpn = self.topydict['DOF_PN'] #  DOF per node
        self.e2sdofmapi = self.topydict['E2SDOFMAPI'] #  Elem to structdof map
        self.nelx = self.topydict['NUM_ELEM_X'] #  Number of elements in X
        self.nely = self.topydict['NUM_ELEM_Y'] #  Number of elements in Y
        self.nelz = self.topydict['NUM_ELEM_Z'] #  Number of elements in Z
        self.fixdof = self.topydict['FIX_DOF'] #  Fixed dof vector
        self.loaddof = self.topydict['LOAD_DOF'] #  Loaded dof vector
        self.loadval = self.topydict['LOAD_VAL'] #  Loaded dof values
        self.Ke = self.topydict['ELEM_K'] #  Element stiffness matrix
        self.Be = self.topydict['ELEM_B'] #  Element strain interpolation matrix
        self.Ce = self.topydict['ELEM_C'] #  Element linear elastic modulus matrix
        self.K = self.topydict['K'] #  Global stiffness matrix
        if self.nelz:
            logger.info('Domain discretisation (NUM_ELEM_X x NUM_ELEM_Y x ' + \
                'NUM_ELEM_Z) = %d x %d x %d' % (self.nelx, self.nely, self.nelz))
        else:
            logger.info( 'Domain discretisation (NUM_ELEM_X x NUM_ELEM_Y) = %d x %d'\
                % (self.nelx, self.nely))

        logger.info('Element type (ELEM_K) = {}'.format(self.topydict['ELEM_TYPE']))
        logger.info('Filter radius (FILT_RAD) = {}'.format(self.filtrad))

        # Check for either one of the following two, will take NUM_ITER if both
        # are specified.
        try:
            self.numiter = self.topydict['NUM_ITER'] #  Number of iterations
            logger.info('Number of iterations (NUM_ITER) = %d' % (self.numiter))
        except KeyError:
            self.chgstop = self.topydict['CHG_STOP'] #  Change stop criteria
            logger.info('Change stop value (CHG_STOP) = %.3e (%.2f%%)' \
                % (self.chgstop, self.chgstop * 100))
            self.numiter = MAX_ITERS

        self.stress = 0.0 # Current maximum stress
        self.objfval = 0.0 # Objective function

        # For stress calculations
        if self.nelz == 0:
            from sympy import lambdify
            from sympy.abc import x, y
            self.Bf = lambdify((x, y), self.Be.copy(), "numpy")
        else:
            from sympy import lambdify
            from sympy.abc import x, y, z
            self.Bf = lambdify((x, y, z), self.Be.copy(), "numpy")

        # All DOF vector and design variables arrays:
        # This needs to be recoded at some point, perhaps. I originally
        # (incorrectly) thought I could do everything just by looking at DOF
        # per node, not so, Cartesian dimension also plays a role.
        # Thus, the 'if'* below is a hack for this to work, and it does...
        if self.dofpn == 1:
            if self.nelz == 0: #  *had to this
                self.e2sdofmapi = self.e2sdofmapi[0:4]
                self.alldof = np.arange(self.dofpn * (self.nelx + 1) * \
                    (self.nely + 1))
                self.desvars = self.volfrac * np.ones((self.nely, self.nelx))
                self.stress_mat = np.zeros((self.nely, self.nelx))
            else:
                self.alldof = np.arange(self.dofpn * (self.nelx + 1) * \
                    (self.nely + 1) * (self.nelz + 1))
                self.desvars = self.volfrac * np.ones((self.nelz, self.nely, self.nelx))
                self.stress_mat = np.zeros((self.nelz, self.nely, self.nelx))
        elif self.dofpn == 2:
            self.alldof = np.arange(self.dofpn * (self.nelx + 1) * (self.nely + 1))
            self.desvars = self.volfrac * np.ones((self.nely, self.nelx))
            self.stress_mat = np.zeros((self.nely, self.nelx))
        else:
            self.alldof = np.arange(self.dofpn * (self.nelx + 1) *\
                (self.nely + 1) * (self.nelz + 1))
            self.desvars = self.volfrac * np.ones((self.nelz, self.nely, self.nelx))
            self.stress_mat = np.zeros((self.nelz, self.nely, self.nelx))
        self.df = np.zeros_like(self.desvars) #  Derivatives of obj. func. (array)
        self.freedof = np.setdiff1d(self.alldof, self.fixdof) #  Free DOF vector
        self.r = np.zeros_like(self.alldof, dtype=np.double) #  Load vector
        self.r[self.loaddof] = self.loadval #  Assign load values at loaded dof
        self.rfree = np.asarray(self.r[self.freedof], dtype=np.double) #  Modified load vector (free dof)
        self.d = np.zeros_like(self.r) #  Displacement vector
        self.dfree = np.zeros_like(self.rfree) #  Modified load vector (free dof)
        # Determine which rows and columns must be deleted from global K:
        self._rcfixed = np.where(np.in1d(self.alldof, self.fixdof), False, True)

        # Print this to screen, just so that the user knows what type of
        # problem is being solved:
        logger.info('Problem type (PROB_TYPE) = ' + self.probtype)
        logger.info('Problem name (PROB_NAME) = ' + self.probname)

        # Set extra parameters if specified:
        # (1) Continuation parameters for 'p':
        self._pmax = self.topydict.get('P_MAX', 1)
        self._phold = self.topydict.get('P_HOLD', self.numiter)
        self._pincr = self.topydict.get('P_INCR')
        self._pcon = self.topydict.get('P_CON', self.numiter)

        # (2) Extra penalisation factor (q) and continuation parameters:
        self.q = self.topydict.get('Q_FAC', 1)

        self._qmax = self.topydict.get('Q_MAX', self.q)
        self._qhold = self.topydict.get('Q_HOLD', self.numiter)
        self._qincr = self.topydict.get('Q_INCR')
        self._qcon = self.topydict.get('Q_CON', self.numiter)


        # (3) Exponential approximation of eta:
        if self.topydict['ETA'] == 'exp':
            #  Initial value of exponent for comp and heat problems:
            self.a = - np.ones(self.desvars.shape)
            if self.probtype == 'mech':
                #  Initial value of exponent for mech problems:
                self.a = self.a * 7 / 3
            self.eta = 1 / (1 - self.a)
            logger.info('Damping factor (ETA) = exp')
        else:
            self.eta = float(self.topydict['ETA']) * np.ones(self.desvars.shape)
            logger.info('Damping factor (ETA) = %3.2f' % (self.eta.mean()))

        try:
            self.approx = self.topydict['APPROX'].lower()
        except KeyError:
            self.approx = None
        if self.approx == 'dquad':
            logger.info('Using diagonal quadratic approximation (APPROX = dquad)')
        # (5) Set passive elements:
        self.pasv = self.topydict['PASV_ELEM']
        if self.pasv.any():
            logger.info('Passive elements (PASV_ELEM) specified')
        else:
            logger.info('No passive elements (PASV_ELEM) specified')
        # (6) Set active elements:
        self.actv = self.topydict['ACTV_ELEM']
        if self.actv.any():
            logger.info('Active elements (ACTV_ELEM) specified')
        else:
            logger.info('No active elements (ACTV_ELEM) specified')

        # Set parameters for compliant mechanism synthesis, if they exist:
        if self.probtype == 'mech':
            if self.topydict['LOAD_DOF_OUT'].any() and \
            self.topydict['LOAD_VAL_OUT'].any():
                self.loaddofout = self.topydict['LOAD_DOF_OUT']
                self.loadvalout = self.topydict['LOAD_VAL_OUT']
            else:
                raise Exception('Not enough input data for mechanism \
                    synthesis!')

            self.rout = np.zeros_like(self.alldof).astype(float)
            self.rout[self.loaddofout] = self.loadvalout
            self.rfreeout = np.array(self.rout[self.freedof], dtype=np.double)
            self.dout = np.zeros_like(self.rout)
            self.dfreeout = np.zeros_like(self.rfreeout)
            ksin = np.ones(self.loaddof.shape, dtype='int') * KDATUM
            ksout = np.ones(self.loaddofout.shape, dtype='int') * KDATUM
            if len(ksin) > 1:
                self.K = self._update_add_mask_sym(self.K, np.asarray([ksin, ksin]), self.loaddof)
                self.K = self._update_add_mask_sym(self.K, np.asarray([ksout, ksout]), self.loaddofout)
            else:
                self.K = self._update_add_mask_sym(self.K, np.asarray([ksin]), self.loaddof)
                self.K = self._update_add_mask_sym(self.K, np.asarray([ksout]), self.loaddofout)

    def fea(self):
        """
        Performs a Finite Element Analysis given the updated global stiffness
        matrix [K] and the load vector {r}, both of which must be in the
        modified state, i.e., [K] and {r} must represent the unconstrained
        system of equations. Return the global displacement vector {d} as a
        NumPy array.

        EXAMPLES:
            >>> t.fea()

        See also: set_top_params

        """
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        if rank == 0:
            if not self.topydict:
                raise Exception('You must first load a TPD file!')
            if self.itercount >= MAX_ITERS:
                raise Exception('Maximum internal number of iterations exceeded!')

            Kfree = self._updateK(self.K.copy())

        # MUMPS with PyMumps
        ctx = DMumpsContext(par=1, sym=1, comm=None)
        ctx.set_icntl(6, 7)
        ctx.set_silent()
        if ctx.myid == 0:
            ctx.set_centralized_sparse(Kfree.tocoo())

        ctx.run(job=4) # Analysis + Factorization

        if ctx.myid == 0:
            x = self.rfree.copy()
            ctx.set_rhs(x) # Modified in place

        ctx.run(job=3) # Solve

        if ctx.myid == 0:
            self.dfree = x

        if self.probtype == 'mech':
            if ctx.myid == 0:
                x = self.rfreeout.copy()
                ctx.set_rhs(x) # Modified in place

            ctx.run(job=3) # Solve

            if ctx.myid == 0:
                self.dfreeout = x

        ctx.destroy() # Cleanup

        # Increment internal iteration counter
        self.itercount += 1

        if rank == 0:
            # Update displacement vectors:
            self.d[self.freedof] = self.dfree
            if self.probtype == 'mech':  # 'adjoint' vectors
                self.dout[self.freedof] = self.dfreeout

            # The following part does not really make sense for heat problems
            if self.probtype == 'heat':
                return

            _L = self.topydict['ELEM_L']
            # Calculate strain and stress values:
            if self.dofpn < 3 and self.nelz == 0:
                Bl = [self.Bf(-_L, -_L), self.Bf(_L, -_L), self.Bf(-_L, _L), self.Bf(_L, _L)]
                for _y in range(self.nely):
                    for _x in range(self.nelx):
                        if self.desvars[_y, _x] > 0:
                            e2sdofmap = self.e2sdofmapi + self.dofpn *\
                                        (_y + _x * (self.nely + 1))
                            stressveclist = []
                            for B in Bl:
                                if self.probtype == 'comp':
                                    strain_vec = np.dot(B, self.d[e2sdofmap])
                                elif self.probtype == 'mech':
                                    strain_vec = np.dot(B, self.d[e2sdofmap] + self.dout[e2sdofmap])

                                stressveclist.append(np.abs((self.desvars[_y, _x] ** self.p) * np.dot(self.Ce, strain_vec)))
                            stress_vec = np.amax(stressveclist, axis=0)
                            self.stress_mat[_y, _x] = self._von_Mises(stress_vec[0], stress_vec[1], 0, stress_vec[2], 0, 0)

            else:
                Bl = [self.Bf(-_L, -_L, -_L), self.Bf(-_L, _L, -_L), self.Bf(-_L, -_L, _L), self.Bf(-_L, _L, _L),
                      self.Bf( _L, -_L, -_L), self.Bf( _L, _L, -_L), self.Bf( _L, -_L, _L), self.Bf( _L, _L, _L)]
                for _z in range(self.nelz):
                    for _y in range(self.nely):
                        for _x in range(self.nelx):
                            if self.desvars[_z, _y, _x] > 0:
                                e2sdofmap = self.e2sdofmapi + self.dofpn *\
                                            (_y + _x * (self.nely + 1) + _z *\
                                            (self.nelx + 1) * (self.nely + 1))
                                stressveclist = []
                                for B in Bl: 
                                    if self.probtype == 'comp':
                                        strain_vec = np.dot(B, self.d[e2sdofmap])
                                    elif self.probtype == 'mech':
                                        strain_vec = np.dot(B, self.d[e2sdofmap] + self.dout[e2sdofmap])

                                    stressveclist.append(np.abs((self.desvars[_z, _y, _x] ** self.p) * np.dot(self.Ce, strain_vec)))
                                stress_vec = np.amax(stressveclist, axis=0)
                                self.stress_mat[_z, _y, _x] = self._von_Mises(stress_vec[0], stress_vec[1], stress_vec[2], stress_vec[3], stress_vec[4], stress_vec[5])

            # Maximum stress reached
            self.stress = np.max(self.stress_mat)


            # Display stress map
            stress_img = self.stress_mat.copy()/self.stress
            if self.dofpn < 3 and self.nelz == 0:
                params = {
                    'prefix': self.probname+"_stress",
                    'iternum': self.itercount,
                    'time': 'none',
                    'filetype': 'png',
                    'dir': "iterations"
                }
                create_2d_imag(stress_img, **params)
            else:
                params = {
                    'prefix': self.probname+"_stress",
                    'iternum': self.itercount,
                    'time': 'none',
                    'dir': "iterations"
                }
                create_3d_geom(stress_img, **params)


    def filter_sens_sigmund(self):
        """
        Filter the design sensitivities using Sigmund's heuristic approach.
        Return the filtered sensitivities.

        EXAMPLES:
            >>> t.filter_sens_sigmund()

        See also: sens_analysis

        """
        if not self.topydict:
            raise Exception('You must first load a TPD file!')
        tmp = np.zeros_like(self.df)
        rmin = int(np.floor(self.filtrad))
        if self.nelz == 0:
            U, V = np.indices((self.nelx, self.nely))
            for i in range(self.nelx):
                umin = np.maximum(i - rmin - 1, 0)
                umax = np.minimum(i + rmin + 2, self.nelx + 1)
                for j in range(self.nely):
                    vmin = np.maximum(j - rmin - 1, 0)
                    vmax = np.minimum(j + rmin + 2, self.nely + 1)
                    u = U[umin: umax, vmin: vmax]
                    v = V[umin: umax, vmin: vmax]
                    dist = self.filtrad - np.sqrt((i - u) ** 2 + (j - v) ** 2)
                    sumnumr = (np.maximum(0, dist) * self.desvars[v, u] *\
                               self.df[v, u]).sum()
                    sumconv = np.maximum(0, dist).sum()
                    tmp[j, i] = sumnumr / (sumconv * self.desvars[j, i])
        else:
            rmin3 = rmin
            U, V, W = np.indices((self.nelx, self.nely, self.nelz))
            for i in range(self.nelx):
                umin, umax = np.maximum(i - rmin - 1, 0),\
                             np.minimum(i + rmin + 2, self.nelx + 1)
                for j in range(self.nely):
                    vmin, vmax = np.maximum(j - rmin - 1, 0),\
                                 np.minimum(j + rmin + 2, self.nely + 1)
                    for k in range(self.nelz):
                        wmin, wmax = np.maximum(k - rmin3 - 1, 0),\
                                     np.minimum(k + rmin3 + 2, self.nelz + 1)
                        u = U[umin:umax, vmin:vmax, wmin:wmax]
                        v = V[umin:umax, vmin:vmax, wmin:wmax]
                        w = W[umin:umax, vmin:vmax, wmin:wmax]
                        dist = self.filtrad - np.sqrt((i - u) ** 2 + (j - v) **\
                               2 + (k - w) ** 2)
                        sumnumr = (np.maximum(0, dist) * self.desvars[w, v, u] *\
                                  self.df[w, v, u]).sum()
                        sumconv = np.maximum(0, dist).sum()
                        tmp[k, j, i] = sumnumr/(sumconv *\
                        self.desvars[k, j, i])

        self.df = tmp

    def sens_analysis(self):
        """
        Determine the objective function value and perform sensitivity analysis
        (find the derivatives of objective function). 

        EXAMPLES:
            >>> t.sens_analysis()

        See also: fea

        """
        if not self.topydict:
            raise Exception('You must first load a TPD file!')
        tmp = self.df.copy()
        self.objfval  = 0.0 #  Objective function value
        
        # Prepare Supporting Variables
        if self.nelz == 0: #  2D problem
    
            Y, X = np.indices((self.nely, self.nelx))
            e2sdofmap = np.expand_dims(self.e2sdofmapi.reshape(-1,1), axis=1)
            e2sdofmap = ne.evaluate("e2sdofmap + dofpn * (Y + X * (nely + 1))", global_dict={"dofpn":self.dofpn, "nely":self.nely})
            Qe = self.d[e2sdofmap]
            QeK = np.tensordot(Qe, self.Ke, axes=(0,0))
            Qe_T = np.swapaxes(Qe, 2, 1).T
            QeKQe = np.einsum('mnk,mnk->mn', QeK, Qe_T)
            
        else: #  3D problem
            
            Z, Y, X = np.indices((self.nelz, self.nely, self.nelx))
            X *= (self.nely + 1)
            Z *= (self.nelx + 1) * (self.nely + 1)
            e2sdofmap = np.expand_dims(self.e2sdofmapi.reshape(-1, 1, 1), axis=1)
            e2sdofmap = ne.evaluate("e2sdofmap + dofpn * (X + Y + Z)", global_dict={"dofpn":self.dofpn})
            Qe = self.d[e2sdofmap]
            QeK = np.tensordot(Qe, self.Ke, axes=(0,0))
            Qe_T = np.swapaxes(Qe.T, 2, 0)
            QeKQe = np.einsum('klmn,klmn->klm', QeK, Qe_T)
        
        # Update TMP
        if self.probtype == 'comp':
            self.objfval += ne.evaluate("sum((desvars ** p) * QeKQe)", global_dict={"desvars":self.desvars, "p":self.p})
            tmp = ne.evaluate("- p * desvars ** (p - 1) * QeKQe", global_dict={"desvars":self.desvars, "p":self.p})

        elif self.probtype == 'heat':
            obj = ne.evaluate("VOID + (1 - VOID) * desvars ** p", global_dict={"desvars":self.desvars, "p":self.p, "VOID":VOID})
            self.objfval += ne.evaluate("sum(obj * QeKQe)")
            fac1 = ne.evaluate("- (1 - VOID) * p * desvars ** (p - 1)", global_dict={"desvars":self.desvars, "p":self.p, "VOID":VOID})
            tmp = ne.evaluate("fac1 * QeKQe")

        elif self.probtype == 'mech':
            self.objfval = ne.evaluate("sum(d)", global_dict={"d":self.d[self.loaddofout]})
            tmp = ne.evaluate("p * desvars ** (p - 1)", global_dict={"desvars":self.desvars, "p":self.p})
            
            if self.nelz == 0:
                QeOut_T = np.swapaxes(self.dout[e2sdofmap], 2, 1).T
                op = np.einsum('mnk,mnk->mn', QeK, QeOut_T)
            else:
                QeOut_T = np.swapaxes(self.dout[e2sdofmap].T, 2, 0)
                op = np.einsum('klmn,klmn->klm', QeK, QeOut_T)
            tmp *= op
            
            
        self.df = tmp


    def update_desvars_oc(self):
        """
        Update the design variables by means of OC-like or equivalently SAO
        method, using the filtered sensitivities; return the updated design
        variables.

        EXAMPLES:
            >>> t.update_desvars_oc()

        See also: sens_analysis, filter_sens_sigmund

        """
        if not self.topydict:
            raise Exception('You must first load a TPD file!')
        # 'p' stays constant for a specified number of iterations from start.
        # 'p' is incremented, but not more than the maximum allowable value.
        # If continuation parameters are not specified in the input file, 'p'
        # will stay constant:
        if self.pcount >= self._phold:
            if (self.p + self._pincr) < self._pmax:
                if (self.pcount - self._phold) % self._pcon == 0:
                    self.p += self._pincr

        if self.qcount >= self._qhold:
            if (self.q + self._qincr) < self._qmax:
                if (self.qcount - self._qhold) % self._qcon == 0:
                    self.q += self._qincr

        self.pcount += 1
        self.qcount += 1

        # Exponential approximation of eta (damping factor):
        if self.itercount > 1:
            if self.topydict['ETA'] == 'exp': #  Check TPD specified value
                mask = np.equal(self.desvarsold / self.desvars, 1)
                self.a = 1 + np.log2(np.abs(self.dfold / self.df)) / \
                np.log2(self.desvarsold / self.desvars + mask) + \
                mask * (self.a - 1)
                self.a = np.clip(self.a, A_LOW, A_UPP)
                self.eta = 1 / (1 - self.a)

        self.dfold = self.df.copy()
        self.desvarsold = self.desvars.copy()

        # Change move limit for compliant mechanism synthesis:
        if self.probtype == 'mech':
            move = 0.1
        else:
            move = 0.2

        lam1, lam2 = 0, 100e3
        dims = self.desvars.shape
        while (lam2 - lam1) / (lam2 + lam1) > 1e-8 and lam2 > 1e-40:
            lammid = 0.5 * (lam1 + lam2)
            if self.probtype == 'mech':
                if self.approx == 'dquad':
                    curv = - 1 / (self.eta * self.desvars) * self.df
                    beta = np.maximum(self.desvars-(self.df + lammid)/curv, VOID)
                    move_upper = np.minimum(move, self.desvars / 3)
                    desvars = np.maximum(VOID, np.maximum((self.desvars - move),\
                    np.minimum(SOLID,  np.minimum((self.desvars + move), \
                    (self.desvars * np.maximum(1e-10, \
                    (-self.df / lammid))**self.eta)**self.q))))
                else:  # reciprocal or exponential
                    desvars = np.maximum(VOID, np.maximum((self.desvars - move),\
                    np.minimum(SOLID,  np.minimum((self.desvars + move), \
                    (self.desvars * np.maximum(1e-10, \
                    (-self.df / lammid))**self.eta)**self.q))))
            else:  # compliance or heat
                if self.approx == 'dquad':
                    curv = - 1 / (self.eta * self.desvars) * self.df
                    beta = np.maximum(self.desvars-(self.df + lammid)/curv, VOID)
                    move_upper = np.minimum(move, self.desvars / 3)
                    desvars = np.maximum(VOID, np.maximum((self.desvars - move),\
                    np.minimum(SOLID,  np.minimum((self.desvars + move_upper), \
                    beta**self.q))))
                else:  # reciprocal or exponential
                    desvars = np.maximum(VOID, np.maximum((self.desvars - move),\
                    np.minimum(SOLID,  np.minimum((self.desvars + move), \
                    (self.desvars * (-self.df / lammid)**self.eta)**self.q))))

            # Check for passive and active elements, modify updated x:
            if self.pasv.any() or self.actv.any():
                flatx = desvars.flatten()
                idx = []
                if self.nelz == 0:
                    y, x = dims
                    for j in range(x):
                        for k in range(y):
                            idx.append(k*x + j)
                else:
                    z, y, x = dims
                    for i in range(z):
                        for j in range(x):
                            for k in range(y):
                                idx.append(k*x + j + i*x*y)
                if self.pasv.any():
                    pasv = np.take(idx, self.pasv) #  new indices
                    np.put(flatx, pasv, VOID) #  = zero density
                if self.actv.any():
                    actv = np.take(idx, self.actv) #  new indices
                    np.put(flatx, actv, SOLID) #  = solid
                desvars = flatx.reshape(dims)

            if self.nelz == 0:
                if desvars.sum() - self.nelx * self.nely * self.volfrac > 0:
                    lam1 = lammid
                else:
                    lam2 = lammid
            else:
                if desvars.sum() - self.nelx * self.nely * self.nelz *\
                self.volfrac > 0:
                    lam1 = lammid
                else:
                    lam2 = lammid
        self.lam = lammid

        self.desvars = desvars

        # Change in design variables:
        self.change = (np.abs(self.desvars - self.desvarsold)).max()

        # Solid-void fraction:
        nr_s = self.desvars.flatten().tolist().count(SOLID)
        nr_v = self.desvars.flatten().tolist().count(VOID)
        self.svtfrac = (nr_s + nr_v) / self.desvars.size

    # ===================================
    # === Private methods and helpers ===
    # ===================================
    def _updateK(self, K):
        """
        Update the global stiffness matrix by looking at each element's
        contribution i.t.o. design domain density and the penalisation factor.
        Return unconstrained stiffness matrix.

        """
        if self.nelz == 0: #  2D problem
            for elx in range(self.nelx):
                for ely in range(self.nely):
                    e2sdofmap = self.e2sdofmapi + self.dofpn *\
                    (ely + elx * (self.nely + 1))
                    if self.probtype == 'comp' or self.probtype == 'mech':
                        updatedKe = self.desvars[ely, elx] ** self.p * self.Ke
                    elif self.probtype == 'heat':
                        updatedKe = (VOID + (1 - VOID) * \
                        self.desvars[ely, elx] ** self.p) * self.Ke
                    K = self._update_add_mask_sym(K, updatedKe, e2sdofmap)
        else: #  3D problem
            for elz in range(self.nelz):
                for elx in range(self.nelx):
                    for ely in range(self.nely):
                        e2sdofmap = self.e2sdofmapi + self.dofpn *\
                                    (ely + elx * (self.nely + 1) + elz *\
                                    (self.nelx + 1) * (self.nely + 1))
                        if self.probtype == 'comp' or self.probtype == 'mech':
                            updatedKe = self.desvars[elz, ely, elx] ** \
                            self.p * self.Ke
                        elif self.probtype == 'heat':
                            updatedKe = (VOID + (1 - VOID) * \
                            self.desvars[elz, ely, elx] ** self.p) * self.Ke
                        K = self._update_add_mask_sym(K, updatedKe, e2sdofmap)


        # Del constrained rows and columns
        # It is faster to convert and operate on the matrix than to slice the
        # original dok_matrix.
        K = K.tocsr()
        K = K[self._rcfixed]
        K = K.tocsc()
        K = K[:,self._rcfixed]

        return K

    # Adapted from the PySparse documentation, in order to act as a substitute
    # for the method of the same name it provided for sparse matrices.
    # Fastest way I could find to do this, but currently it is the
    # implementation's greatest bottleneck by far.
    # Fills just the upper triangle of the matrix, because it is symmetric and
    # it fits with MUMPS's way of handling sym-pos-def matrices.
    def _update_add_mask_sym(self, A, B, ind):
        for i in range(len(ind)):
            for j in range(i, len(ind)):
                k = dict.get(A, (ind[i], ind[j]), 0)
                A._update({(ind[i],ind[j]):k+B[i,j]})

        return A

    # Calculates von Mises criterion
    def _von_Mises(self, s11, s22, s33, s12, s13, s23):
        return np.sqrt(0.5*((s11 - s22)**2 + (s22 - s33)**2 + (s33 - s11)**2 + 6*(s12**2 + s13**2 + s23**2)))


# EOF topology.py
