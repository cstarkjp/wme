"""
---------------------------------------------------------------------

Module applying the weathering-mediated erosion model to bedrock channel geometry

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
    Numerical solution of combined models of weathering-mediated erosion & bedrock channel
    
    Class that provides numerical solution of a combination model 1d weathering-mediated
    erosion and 1+1d bedrock channel cross-section (at the channel wall).
    
    Args:
        em (:class:`~.theory.WeatheringMediatedErosion`): 
                instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        pdict (dict): model parameters dictionary
        
        
    """
    def __init__(self, em, pdict):
        """
        Initialize class instance.

        Attributes:
            pdict (:obj:`dict`) : 
                model parameters dictionary, extended during & after instantiation

            v0_eqn_vr_h_z (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                
            w0_eqn_wr_z (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                
            vs_eqn_vr (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                
            w0_calibrated (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                
            v0_calibrated (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                
            W_calibrated  (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                
            vs_calibrated (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                

            
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
        Compute dependence of various properties with height above channel base

        Args:
            n_pts (int): 
                number of sampling points along vertical profile  

        Attributes:
            n_pts (:obj:`int`): 
                record of number of sampling points
            z_array  (:class:`numpy.ndarray`) :
                sample points along vertical profile 
            vs_calibrated (:class:`sympy.Eq <sympy.core.relational.Equality>`) : 
                weathering-mediated surface-normal erosion rate at steady state
                calibrated by model parameters supplied in :attr:`wme.pdict`
            w0_array (:class:`numpy.ndarray`) :
                surface weakness along vertical profile (at steady state)
            v0_array (:class:`numpy.ndarray`) :
                baseline (absent weathering-driven weakening) surface-normal erosion rate 
                along vertical profile  (at steady state)
            vs_array (:class:`numpy.ndarray`) :
                weathering-mediated surface-normal erosion rate along vertical profile 
                (at steady state)
            eta0_array (:class:`numpy.ndarray`) :
                surface weakness along vertical profile (at steady state)
            W_array (:class:`numpy.ndarray`) :
                weathering number along vertical profile

        """
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
        Compute dependence of various properties with height above channel base
        
        Attributes:
            vs_array (:class:`numpy.ndarray`) :
                
            dzdy_array (:class:`numpy.ndarray`) :
                
            y_array (:class:`numpy.ndarray`) :
                
            ch_y_array (:class:`numpy.ndarray`) :
                
            ch_z_array (:class:`numpy.ndarray`) :
                

        """
        self.vs_array = np.array([
            np.float64(self.vs_calibrated.rhs.subs(z,z__)) for z__ in self.z_array])
        
        self.dzdy_array = 1/np.sqrt(((np.max(self.vs_array)+0.01)/self.vs_array)**2-1)
        self.y_array = integrate.cumtrapz(self.dzdy_array,x=self.z_array, initial=0)
        self.ch_y_array = np.concatenate(
            (np.linspace(-np.max(self.z_array)*0.3,0,self.z_array.shape[0]),self.y_array))
        self.ch_z_array \
            = np.concatenate(( np.zeros(self.y_array.shape[0]),self.z_array ))
    
        
        
        
        