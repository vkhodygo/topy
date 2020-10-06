"""
# =============================================================================
# A class to optimise the topology of a design domain for defined boundary
# conditions. Data is read from an input file, see 'examples' directory.
#
# Authors: William Hunter, Tarcísio L. de Oliveira
# Copyright (C) 2008, 2015, William Hunter.
# Copyright (C) 2020, Tarcísio L. de Oliveira
# =============================================================================
"""
import os

import numpy as np
#from pysparse import superlu, itsolvers, precon
from scipy import sparse
from scipy.sparse import linalg
from sympy import symbols

from .utils import get_logger
from .parser import tpd_file2dict, config2dict
from .visualisation import *
from .pathfinding import *

logger = get_logger(__name__)
__all__ = ['TopologyGen']


MAX_ITERS = 250

SOLID, VOID = 1.000, 0.001 #  Upper and lower bound value for design variables
KDATUM = 0.1 #  Reference stiffness value of springs for mechanism synthesis

# Constants for exponential approximation:
A_LOW = -3 #  Lower restriction on 'a' for exponential approximation
A_UPP = -1e-5 #  Upper restriction on 'a' for exponential approximation


# =======================
# === ToPy base class ===
# =======================
class TopologyGen:
    """
    A class to optimise the topology of a design domain for defined boundary
    values. Data is read from an input file (see 'examples' folder).

    """
    def __init__(self, config=None, topydict={}, pcount=0, 
                 qcount=0, itercount=0, change=1, svtfrac=None):
        self.pcount = pcount #  Counter for continuation of p
        self.qcount = qcount #  Counter for continuation of q for GSF
        self.itercount = itercount-1 #  Internal counter
        self.change = change
        self.svtfrac = svtfrac

        if config:
            self.topydict = config2dict(config.copy())
        else:
            self.topydict = topydict #  Store tpd file data in dictionary


    # ======================
    # === Public methods ===
    # ======================
    def preprocess_space(self):
        Kfree = self.createK()
        self.fea(Kfree)
        if self.probtype != "heat":
            logger.info("\nBase stress: %3.1f\n" % (self.stress*1e-6))
            self.expand()
        self.desvars = self.volfrac*self.desvars
        Kfree = self.preprocessK()

        return Kfree


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
        self.Smax = self.topydict['STRESS_MAX'] # Maximum stress allowed
        if self.nelz:
            logger.info('Domain discretisation (NUM_ELEM_X x NUM_ELEM_Y x ' + \
                'NUM_ELEM_Z) = %d x %d x %d' % (self.nelx, self.nely, self.nelz))
        else:
            logger.info( 'Domain discretisation (NUM_ELEM_X x NUM_ELEM_Y) = %d x %d'\
                % (self.nelx, self.nely))

        if type(self.loaddof) is not list:
            self.loaddof = np.asarray(self.loaddof)

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
        self.strain = 0.0 # Current maximum strain
        self.disp   = 0.0 # Current maximum displacement

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
                self.desvars = np.zeros((self.nely, self.nelx))
                self.stress_mat = np.zeros((self.nely, self.nelx))
            else:
                self.alldof = np.arange(self.dofpn * (self.nelx + 1) * \
                    (self.nely + 1) * (self.nelz + 1))
                self.desvars = np.zeros((self.nelz, self.nely, self.nelx))
                self.stress_mat = np.zeros((self.nelz, self.nely, self.nelx))
        elif self.dofpn == 2:
            self.alldof = np.arange(self.dofpn * (self.nelx + 1) * (self.nely + 1))
            self.desvars = np.zeros((self.nely, self.nelx))
            self.stress_mat = np.zeros((self.nely, self.nelx))
        else:
            self.alldof = np.arange(self.dofpn * (self.nelx + 1) *\
                (self.nely + 1) * (self.nelz + 1))
            self.desvars = np.zeros((self.nelz, self.nely, self.nelx))
            self.stress_mat = np.zeros((self.nelz, self.nely, self.nelx))
        self.df = np.zeros_like(self.desvars) #  Derivatives of obj. func. (array)
        self.freedof = np.setdiff1d(self.alldof, self.fixdof) #  Free DOF vector
        self.r = np.zeros_like(self.alldof, dtype=np.float32) #  Load vector
        self.r[self.loaddof] = self.loadval #  Assign load values at loaded dof
        self.rfree = np.asarray(self.r[self.freedof], dtype=np.float32) #  Modified load vector (free dof)
        self.d = np.zeros_like(self.r) #  Displacement vector
        self.dfree = np.zeros_like(self.rfree) #  Modified load vector (free dof)
        # Determine which rows and columns must be deleted from global K:
        self._rcfixed = np.where(np.in1d(self.alldof, self.fixdof), False, True)
        self.remove_list = self._rcfixed

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
            self.rfreeout = np.array(self.rout[self.freedof], dtype=np.float32)
            self.dout = np.zeros_like(self.rout)
            self.dfreeout = np.zeros_like(self.rfreeout)
            ksin = np.ones(self.loaddof.shape, dtype='int') * KDATUM
            ksout = np.ones(self.loaddofout.shape, dtype='int') * KDATUM
            maskin = np.ones(self.loaddof.shape, dtype='int')
            maskout = np.ones(self.loaddofout.shape, dtype='int')
            if len(ksin) > 1:
                self.K = self._update_add_mask_sym(self.K, np.asarray([ksin, ksin]), self.loaddof, maskin)
                self.K = self._update_add_mask_sym(self.K, np.asarray([ksout, ksout]), self.loaddofout, maskout)
            else:
                self.K = self._update_add_mask_sym(self.K, np.asarray([ksin]), self.loaddof, maskin)
                self.K = self._update_add_mask_sym(self.K, np.asarray([ksout]), self.loaddofout, maskout)

    def fea(self, Kfree):
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
        if not self.topydict:
            raise Exception('You must first load a TPD file!')
        if self.itercount >= MAX_ITERS:
            raise Exception('Maximum internal number of iterations exceeded!')

        # Kfree = self._updateK(self.K.copy())

        removed = self.remove_list[self._rcfixed].copy()

        if self.dofpn < 3 and self.nelz == 0: #  Direct solver
            Kfree = sparse.csc_matrix(Kfree)
            lu = linalg.splu(Kfree)
            self.dfree[removed] = lu.solve(self.rfree[removed])
            if self.probtype == 'mech':
                self.dfreeout[removed] = lu.solve(self.rfreeout[removed])
        else: #  Iterative solver for 3D problems
            Kfree = sparse.csc_matrix(Kfree, dtype=np.float32)
            lu = linalg.splu(Kfree)
            preK_x = lambda x: lu.solve(np.asarray(x, dtype=np.float32))
            preK = linalg.LinearOperator(Kfree.shape, preK_x)
            # ToPy's original implementation used Preconditioned Conjugate Gradient (PCG) from PySparse.
            # SciPy's implementation (linalg.cg) didn't converge. This was the first function that
            # managed to converge, but it is pretty slow.
            self.dfree, info = linalg.minres(Kfree, self.rfree, tol=1e-8, maxiter=8000, M=preK)
            if info > 0:
                logger.error('Solver error: Number of iterations: {}.'.format(info))
                raise Exception('Solution for FEA did not converge.')
            else:
                logger.debug('TgPy: Solution for FEA converged.')
            if self.probtype == 'mech':  # mechanism synthesis
                # This was also PCG.
                self.dfreeout, info = linalg.minres(Kfree, self.rfreeout, tol=1e-8, maxiter=8000, M=preK)
                if info > 0:
                    logger.error('Solver error: Number of iterations: {}.'.format(info))
                    raise Exception('Solution for FEA of adjoint load '
                                    'case did not converge.')

        # Update displacement vectors:
        self.d[self.freedof] = self.dfree
        if self.probtype == 'mech':  # 'adjoint' vectors
            self.dout[self.freedof] = self.dfreeout

        if self.probtype == 'heat':
            self.itercount += 1
            return

        # Calculate strain and stress values:
        _L = self.topydict['ELEM_L']
        if self.dofpn < 3 and self.nelz == 0:
            x, y = symbols("x y")
            num_elem = self.nelx * self.nely
            for i in range(num_elem):
                _x = int(np.floor(i / self.nely))
                _y = i % self.nely
                e2sdofmap = self.e2sdofmapi + self.dofpn *\
                            (_y + _x * (self.nely + 1))
                if self.desvars[_y, _x] > 0:
                    B = np.array(self.Be.copy().subs({x:(2*_L*(_x+0.5)), y:(2*_L*(_y+0.5))})).astype('double')
                    if self.probtype == 'comp' or self.probtype == 'heat':
                        strain_vec = np.dot(B, self.d[e2sdofmap])
                    elif self.probtype == 'mech':
                        strain_vec = np.dot(B, self.d[e2sdofmap])

                    stress_vec = self.desvars[_y, _x] * np.dot(self.Ce, strain_vec)
                    if self.probtype == 'comp' or self.probtype == 'mech':
                        self.stress_mat[_y, _x] = self._von_Mises(stress_vec[0], stress_vec[1], 0, stress_vec[2], 0, 0)
                    elif self.probtype == 'heat':
                        self.stress_mat[_y, _x] = np.linalg.norm(stress_vec)

        else:
            x, y, z = symbols("x y z")
            num_elem = self.nelx * self.nely
            for i in range(num_elem):
                _z = int(np.floor(i / (self.nely * self.nelx)))
                rest = i % (self.nely * self.nelx)
                _x = int(np.floor(rest / self.nely))
                _y = rest % self.nely
                e2sdofmap = self.e2sdofmapi + self.dofpn *\
                            (ely + elx * (self.nely + 1) + elz *\
                            (self.nelx + 1) * (self.nely + 1))
                if self.desvars[_z, _y, _x] > 0:
                    B = np.array(self.Be.copy().subs({x:(2*_L*(_x+0.5)), y:(2*_L*(_y+0.5)), z:(2*_L*(_z+0.5)) })).astype('double')
                    strain_vec = np.dot(B, self.d[e2sdofmap])
                    stress_vec = self.desvars[_z, _y, _x] * np.dot(self.Ce, strain_vec)
                    self.stress_mat[_z, _y, _x] = self._von_Mises(stress_vec[0], stress_vec[1], stress_vec[2], stress_vec[3], stress_vec[4], stress_vec[5])

        # # Filter force singularities
        # if self.itercount >= 0:
        #     self._saint_venant()
        
        # Maximum stress reached
        self.stress = np.max(self.stress_mat)

        # Increment internal iteration counter
        self.itercount += 1

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
                    if self.desvars[j, i] > 0:
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
                        if self.desvars[k, j, i] > 0:
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
            e2sdofmap = np.add(e2sdofmap, self.dofpn * (Y + X * (self.nely + 1)))
            Qe = self.d[e2sdofmap]
            QeK = np.tensordot(Qe, self.Ke, axes=(0,0))
            Qe_T = np.swapaxes(Qe, 2, 1).T
            QeKQe = np.einsum('mnk,mnk->mn', QeK, Qe_T)
            
        else: #  3D problem
            
            Z, Y, X = np.indices((self.nelz, self.nely, self.nelx))
            X *= (self.nely + 1)
            Z *= (self.nelx + 1) * (self.nely + 1)
            e2sdofmap = np.expand_dims(self.e2sdofmapi.reshape(-1, 1, 1), axis=1)
            e2sdofmap = np.add(e2sdofmap, self.dofpn * (X + Y + Z))
            Qe = self.d[e2sdofmap]
            QeK = np.tensordot(Qe, self.Ke, axes=(0,0))
            Qe_T = np.swapaxes(Qe.T, 2, 0)
            QeKQe = np.einsum('klmn,klmn->klm', QeK, Qe_T)
        
        # Update TMP
        if self.probtype == 'comp':
            self.objfval += ((self.desvars ** self.p) * QeKQe).sum()
            tmp = - self.p * self.desvars ** (self.p - 1) * QeKQe

        elif self.probtype == 'heat':
            obj = (VOID + (1 - VOID) * self.desvars ** self.p)
            self.objfval += (obj * QeKQe).sum()
            fac1 = - (1 - VOID) * self.p * self.desvars ** (self.p - 1)
            tmp = fac1 * QeKQe

        elif self.probtype == 'mech':
            self.objfval = self.d[self.loaddofout].sum()
            tmp = self.p * self.desvars ** (self.p - 1)
            
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

        removed_mask = np.equal(self.desvars, 0)

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

        desvars[removed_mask] = 0
        self.desvars = desvars

        # Change in design variables:
        self.change = (np.abs(self.desvars - self.desvarsold)).max()

        # Solid-void fraction:
        nr_s = self.desvars.flatten().tolist().count(SOLID)
        nr_v = self.desvars.flatten().tolist().count(VOID)
        nr_0 = self.desvars.flatten().tolist().count(0)
        self.svtfrac = (nr_s + nr_v + nr_0) / self.desvars.size

    # ===================================
    # === Private methods and helpers ===
    # ===================================
    def createK(self):
        """
        Update the global stiffness matrix by looking at each element's
        contribution i.t.o. design domain density and the penalisation factor.
        Return unconstrained stiffness matrix.

        """
        K = self.K.copy()
        self.remove_list = self._rcfixed.copy()
        
        loads = split_loads(self.loaddof, self.nelx, self.nely, self.nelz, self.dofpn)
        fixed = split(self.fixdof, self.nelx, self.nely, self.nelz, self.dofpn)
        active = active_points(self.actv, self.nelx, self.nely, self.nelz)
        passive = passive_points(self.pasv, self.nelx, self.nely, self.nelz)

        for l in loads:
            self.desvars = get_path(self.desvars, fixed, active, passive, l)

        if self.probtype == 'mech':
            loads_out = split_loads(self.loaddofout, self.nelx, self.nely, self.nelz, self.dofpn)
            for lo in loads_out:
                self.desvars = get_path(self.desvars, loads, active, passive, lo)
                self.desvars = get_path(self.desvars, fixed, active, passive, lo)

        if self.nelz == 0: #  2D problem
            for elx in range(self.nelx):
                for ely in range(self.nely):
                    if self.desvars[ely, elx] == 1:
                        e2sdofmap = self.e2sdofmapi + self.dofpn *\
                                    (ely + elx * (self.nely + 1))
                        updatedKe = self.Ke
                        mask = np.ones(e2sdofmap.size, dtype=int)
                        K = self._update_add_mask_sym(K, updatedKe, e2sdofmap, mask)

        else: #  3D problem
            for elz in range(self.nelz):
                for elx in range(self.nelx):
                    for ely in range(self.nely):
                        if self.desvars[elz, ely, elx] == 1:
                            e2sdofmap = self.e2sdofmapi + self.dofpn *\
                                        (ely + elx * (self.nely + 1) + elz *\
                                        (self.nelx + 1) * (self.nely + 1))
                            updatedKe = self.Ke
                            mask = np.ones(e2sdofmap.size, dtype=int)
                            K = self._update_add_mask_sym(K, updatedKe, e2sdofmap, mask)

        for i in range(K.shape[0]):
            if K[i, i] == 0:
                self.remove_list[i] = False

        #  Del constrained rows and columns
        K = K[self.remove_list][:,self.remove_list]
        return K

    def updateK(self):
        """
        Update the global stiffness matrix by looking at each element's
        contribution i.t.o. design domain density and the penalisation factor.
        Return unconstrained stiffness matrix.

        """

        K = self.K.copy()
        if self.nelz == 0: #  2D problem
            for elx in range(self.nelx):
                for ely in range(self.nely):
                    if self.desvars[ely, elx] > 0:
                        e2sdofmap = self.e2sdofmapi + self.dofpn *\
                        (ely + elx * (self.nely + 1))
                        if self.probtype == 'comp' or self.probtype == 'mech':
                            updatedKe = self.desvars[ely, elx] ** self.p * self.Ke
                        elif self.probtype == 'heat':
                            updatedKe = (VOID + (1 - VOID) * \
                            self.desvars[ely, elx] ** self.p) * self.Ke
                        mask = np.ones(e2sdofmap.size, dtype=int)
                        K = self._update_add_mask_sym(K, updatedKe, e2sdofmap, mask)
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
                        mask = np.ones(e2sdofmap.size, dtype=int)
                        K = self._update_add_mask_sym(K, updatedKe, e2sdofmap, mask)


        #  Del constrained rows and columns
        K = K[self.remove_list][:,self.remove_list]
        return K

    def expand(self):
        """
        Generate a basic structure based on the preliminary stress analysis
        """
        warn_width = False
        warn_height = False
        if self.nelz == 0: #  2D problem
            for elx in range(self.nelx):
                for ely in range(self.nely):
                    if self.stress_mat[ely, elx] > 0:
                        h = np.ceil(np.sqrt(self.stress_mat[ely, elx]/self.Smax))
                        h += h % 2 # Make h even
                        umin = int(np.maximum(elx - h/2, 0))
                        umax = int(np.minimum(elx + h/2 + 1, self.nelx))
                        if umax - umin < h and not warn_width:
                            logger.info("Warning: design space may have insufficient width.")
                            warn_width = True
                        vmin = int(np.maximum(ely - h/2, 0))
                        vmax = int(np.minimum(ely + h/2 + 1, self.nely))
                        if vmax - vmin < h and not warn_height:
                            logger.info("Warning: design space may have insufficient height.")
                            warn_height = True
                        for i in range(umin, umax):
                            for j in range(vmin, vmax):
                                if (i-elx)**2 + (j-ely)**2 <= (h/2)**2:
                                    self.desvars[j, i] = SOLID
                        #print(self.desvars[vmin:vmax, umin:umax])
        else: #  3D problem
            for elz in range(self.nelz):
                for elx in range(self.nelx):
                    for ely in range(self.nely):
                        print()

    def preprocessK(self):
        """
        Remove elements which will not be used in the problem.
        """
        K = self.K.copy()
        self.remove_list = self._rcfixed.copy()
        if self.nelz == 0: #  2D problem
            for elx in range(self.nelx):
                for ely in range(self.nely):
                    e2sdofmap = self.e2sdofmapi + self.dofpn *\
                    (ely + elx * (self.nely + 1))
                    if self.desvars[ely, elx] > 0:
                        updatedKe = self.Ke
                        mask = np.ones(e2sdofmap.size, dtype=int)
                        K = self._update_add_mask_sym(K, updatedKe, e2sdofmap, mask)

            for i in range(K.shape[0]):
                if K[i, i] == 0:
                    self.remove_list[i] = False

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
                        mask = np.ones(e2sdofmap.size, dtype=int)
                        K = self._update_add_mask_sym(K, updatedKe, e2sdofmap, mask)

        #  Del constrained rows and columns
        K = K[self.remove_list][:,self.remove_list]
        return K

    def _saint_venant(self):
        multipliers = np.asarray([[1.0, 0.5, 0.5, 1.0],
                                  [0.5, 0.0, 0.0, 0.5],
                                  [0.5, 0.0, 0.0, 0.5],
                                  [1.0, 0.5, 0.5, 1.0]])

        if self.nelz == 0: #  2D problem
            nodes = np.floor(self.loaddof / self.dofpn) + 1
            node_x = (np.floor(nodes / (self.nely + 1))).astype(int)
            node_y = (nodes % (self.nely + 1)).astype(int) - 1

            umin = np.maximum(node_x - 2, 0)
            umax = np.minimum(node_x + 2, self.nelx)
            vmin = np.maximum(node_y - 2, 0)
            vmax = np.minimum(node_y + 2, self.nely)
            true_xmin = node_x - 2
            true_ymin = node_y - 2

            for i in range(len(umin)):
                xmin = umin[i]
                xmax = umax[i]
                ymin = vmin[i]
                ymax = vmax[i]

                xmin_m = xmin - true_xmin[i] 
                xmax_m = xmax - true_xmin[i] 
                ymin_m = ymin - true_ymin[i] 
                ymax_m = ymax - true_ymin[i] 

                subsm = self.stress_mat[ymin:ymax, xmin:xmax]
                submult = multipliers[ymin_m:ymax_m, xmin_m:xmax_m]

                self.stress_mat[ymin:ymax, xmin:xmax] = np.multiply(subsm, submult)

        else: #  3D problem
            print()

    # Taken from the PySparse documentation, in order to act as a substitute
    # for the method of the same name it provided for sparse matrices.
    def _update_add_mask_sym(self, A, B, ind, mask):
        for i in range(len(ind)):
            for j in range(len(ind)):
                if mask[i]:
                    A[ind[i],ind[j]] += B[i,j]

        return A

    # Calculates von Mises criterion
    def _von_Mises(self, s11, s22, s33, s12, s13, s23):
        return np.sqrt(0.5*((s11 - s22)**2 + (s22 - s33)**2 + (s33 - s11)**2 + 6*(s12**2 + s13**2 + s23**2)))

# EOF topology.py
