"""
---------------------------------------------------------------------

Module for ... 

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`numpy`
  -  :mod:`sympy`

Imports symbols from :mod:`.symbols` module

---------------------------------------------------------------------

"""


import numpy as np
import sympy as sy
from sympy import Eq
from .symbols import *



class ErosionWeatheringModel:
    """
    TBD
    
    Attributes:
        xx (xx) : xxxx

    """
    def __init__(self):
        """
        Initialize class instance.
        """
        self.W_eqn           = Eq(W,w_0/(k*v_0))
        self.nus_eqn_W       = Eq(nu_s,(1+sy.sqrt(1+4*W))/2)
        self.nus_eqn_w0_v0   = self.nus_eqn_W.subs(W,self.W_eqn.rhs)
        self.etas0_eqn_W     = Eq(eta_s0, self.nus_eqn_W.rhs)
        self.etas0_eqn_w0_v0 = self.etas0_eqn_W.subs(W,self.W_eqn.rhs)
        self.nus_eqn_etas0   = Eq(nu_s, eta_s0)
        self.vs_eqn_etas0_v0 = Eq(v_s, self.nus_eqn_etas0.rhs*v_0)
        self.v0_eqn_etas0_vs = Eq(v_0, sy.solve(self.vs_eqn_etas0_v0,v_0)[0])
        self.vs_eqn_w0_v0    = Eq(v_s, self.nus_eqn_w0_v0.rhs*v_0)
        self.v0_eqn_vs_w0    = Eq(v_0, sy.solve(self.vs_eqn_w0_v0,v_0)[0])
        self.v0_eqn_vr_h_z   = Eq(v_0, v_r*((h-self.step(z,v_smooth,v_off))*(1-v_b)+v_b))
        self.w0_eqn_wr_z     = Eq(w_0, w_r*self.step(z,w_smooth,w_off))


    @staticmethod
    def step(x,k,x0):
        """
        TBD.
        """
        return (1+sy.tanh((x-x0)*k))/2



        