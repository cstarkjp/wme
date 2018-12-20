"""
---------------------------------------------------------------------

Module for 

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`numpy`
  -  :mod:`scipy.integrate` (integrate)
  
Imports symbols from :mod:`.symbols` module

---------------------------------------------------------------------

"""


import numpy as np
from scipy import integrate
from .symbols import *



class ChannelWall:
    """
    Numerical solution of ...
    
    Class that provides ...
    
    
    Args:
        pdict (dict): model parameters dictionary
        
    Attributes:
        pdict (dict) : model parameters dictionary, extended during & after instantiation
        
        
    """
    def __init__(self, em, pdict):
        """
        Initialize class instance.
        """
        self.pdict = pdict
        self.v0_eqn_vr_h_z = em.v0_eqn_vr_h_z.subs(pdict)
        self.w0_eqn_wr_z = em.w0_eqn_wr_z.subs(pdict)
        self.vs_eqn_vr = em.vs_eqn_w0_v0.subs({v_0:self.v0_eqn_vr_h_z.rhs})
                        
        self.w0_calibrated = self.w0_eqn_wr_z.subs(self.pdict)
        self.v0_calibrated = self.v0_eqn_vr_h_z.subs(self.pdict)
        self.W_calibrated  = em.W_eqn.subs(self.pdict)
        self.vs_calibrated \
            = self.vs_eqn_vr.subs({w_0:self.w0_eqn_wr_z.rhs}).subs(self.pdict)

    def compute_vertical_profiles(self, n_pts=100):
        """
        Use ...
        """
#         self.pdict.update({})
        self.n_pts = n_pts
        self.z_array  = np.linspace(0,1,self.n_pts)
        self.w0_array = np.array([
            np.float64(self.w0_calibrated.rhs.subs(z,z__)) for z__ in self.z_array])
        self.v0_array = np.array([
            np.float64(self.v0_calibrated.rhs.subs(z,z__)) for z__ in self.z_array])
        self.vs_array = np.array([
            np.float64(self.vs_calibrated.rhs.subs(z,z__)) for z__ in self.z_array])
        self.eta0_array = self.vs_array/self.v0_array
        self.W_array = np.array([
            np.float64(  (self.W_calibrated.rhs.subs(z,z__))
                            .subs({v_0:self.v0_array[idx]})
                            .subs({w_0:self.w0_array[idx]})  )
                                for idx,z__ in enumerate(self.z_array)])
        
        self.vs_calibrated \
            = self.vs_eqn_vr.subs({w_0:self.w0_eqn_wr_z.rhs}).subs(self.pdict)
    
    
    def compute_cross_section(self):
        """
        Use ...
        """
        self.vs_array = np.array([
            np.float64(self.vs_calibrated.rhs.subs(z,z__)) for z__ in self.z_array])
        
        self.dzdy_array = 1/np.sqrt(((np.max(self.vs_array)+0.01)/self.vs_array)**2-1)
        self.y_array = integrate.cumtrapz(self.dzdy_array,x=self.z_array, initial=0)
        self.ch_y_array = np.concatenate(
            (np.linspace(-np.max(self.z_array)*0.3,0,self.z_array.shape[0]),self.y_array))
        self.ch_z_array \
            = np.concatenate(( np.zeros(self.y_array.shape[0]),self.z_array ))
    
        
        
        
        