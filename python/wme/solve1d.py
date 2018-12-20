"""
---------------------------------------------------------------------

Module for performing finite-difference solution of weathering-driven
weakening of a 1d bedrock surface half-space and its concomitant erosion.

The moving-boundary problem is solved using a continuous-valued 
distance (from the erosion front) function that tracks motion with sub-grid
resolution. Upwind differencing and explicit Euler methods are employed. 

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`numpy`

Imports symbols from :mod:`.symbols` module

---------------------------------------------------------------------

"""


import numpy as np
from .symbols import *


def negExpH(chi,dchi):
    """
    Negate, exponentiate and Heaviside clip.
    
    Assumes the argument chi= :math:`\\chi`
    is the dimensionless distance from the erosion front, 
    and dchi= :math:`\\Delta\\chi`  is the discrete spatial step size.
    Exponentiates -chi*H(chi+dchi) where H=Heaviside function.
    When invoked in computing d{eta}/d{chi}, adding dchi means a non-clipped
    value is returned for the sample point just to the left (-ve chi).
    Thus the eta gradient is estimated across the erosion front, where sample points
    span the moving origin, as well as for all points chi>=0.
    
    Args:
        chi (float): distance from the erosion front
        dchi (float): discrete spacing between sample points along x
    
    Returns:
        float:
            exp(-chi) if chi+dchi>=0, 1 otherwise
    """
    return np.exp(-chi*np.heaviside(chi+dchi,0))

def nu_s_W(W):
    """
    Dimensionless steady-state speed of the erosion front :math:`\\nu_s(W)`.
    
    Assumes: :math:`\\nu_s = \\tfrac{1}{2}\\left(1+\sqrt{1+4W}\\right)`
    
    Args:
        W (float): weathering number :math:`W`
    
    Returns:
        float: dimensionless erosion rate :math:`\\nu_s`
    """
    return 0.5*(1+np.sqrt(1+4*W))

def eta_chi_tau(chi,tau,W):
    """
    Weathering-driven weakness function.
    
    Analytic solution for eta(chi,tau) as a function of dimensionless distance
    (depth into the rock) chi and time tau, and parameterized
    by weathering number W, assuming an exponential-decay model for weathering.
    
    Args:
        tau (float): dimensionless time
        chi (float): dimensionless distance
        W (float): weathering number
    
    Returns:
        float:
            eta(chi,tau;W)
    """
    return (1+(W/nu_s_W(W))*np.exp(-(chi)))*np.heaviside(chi,0)


class ErosionWeathering:
    """
    Numerical solution of eta(chi,tau)  and phi(tau) evolution 
    
    Class that provides a finite-difference method for solving the (chi,tau) evolution
    of a weakness profile eta(chi,tau) and its eroding surface position phi(t)
    as 2d array eta_i^j and 1d vector phi^j respectively.
    and that provides dictionaries for the model and its numerical solution
    parameters.
    
    
    Args:
        pdict (dict): model parameters dictionary
        ndict (dict): numerical method parameters dictionary
        
    Attributes:
        pdict (dict) : model parameters dictionary, extended during & after instantiation
        ndict (dict) : numerical method parameters dictionary
        
        chi_domain_size (float):  length of chi solution domain 
                                    (extracted from ndict)
        Delta_chi (float):        spacing between discrete chi solution points 
                                    (extracted from ndict)
        n_chi_domain (int):       number of solution points in distance chi 
                                    (extracted from ndict)
        tau_domain_size (float):  maximum duration of solution (truncated if/when front 
                                  exits chi domain) (extracted from ndict)
        tau_n_steps (int):        number of solution points in time tau
        Delta_tau (float):        spacing between discrete tau solution points

        chi_array (numpy.ndarray) : chi_i discrete distances
        tau_array (numpy.ndarray) : tau^j discrete times
        eta_array (numpy.ndarray) : eta_i^j discretized weakness profile
        phi_array (numpy.ndarray) : phi^j discrete (in time) series of erosion front 
                                    positions (smoothly resolved as floats)
        nu_array  (numpy.ndarray) : nu^j  discrete (in time) series of 
                                    dimensionless erosion rates
        
        j (int)        :  final time step index

        W (float)    : weathering number
        nu_s (float) : predicted (by analytical solution) dimensionless 
                            steady-state erosion rate
        v_s (float)  : predicted (by analytical solution)  steady-state erosion rate
        nu_s_bar (float) : post-hoc estimate (from averaging portion of solutions) of
                           dimensionless steady-state erosion rate
    """
    def __init__(self, pdict, ndict):
        """
        Initialize class instance.
        """
        self.pdict = pdict
        self.W     = pdict[w_0]/(pdict[k]*pdict[v_0])
        self.nu_s  = nu_s_W(self.W)
        self.v_s   = self.nu_s*pdict[v_0]
        self.pdict.update({W:self.W, nu_s:self.nu_s, v_s:self.v_s})
        
        self.ndict = ndict
        self.chi_domain_size = ndict[chi_domain_size]
        self.tau_domain_size = ndict[tau_domain_size]
        self.Delta_chi = ndict[Delta_chi]
        self.Delta_tau = ndict[Delta_tau]
        self.n_chi_domain = np.int64(self.chi_domain_size/self.Delta_chi)+1
        self.tau_n_steps  = np.int64(self.tau_domain_size/self.Delta_tau)+1

        self.eta_array = np.zeros((self.tau_n_steps,self.n_chi_domain),dtype=np.float64)
        self.eta_array[0] = np.ones(self.n_chi_domain,dtype=np.float64)
        self.phi_array = np.zeros((self.tau_n_steps),dtype=np.float64)
        self.chi_array = np.linspace(0,self.chi_domain_size,self.n_chi_domain,
                                     dtype=np.float64)
        self.tau_array = np.linspace(0,self.tau_domain_size,self.tau_n_steps,
                                     dtype=np.float64)
        self.nu_array  = np.zeros((self.tau_n_steps),dtype=np.float64)
        self.j         = 0
        

    def solve(self):
        """
        Use an explicit finite-difference scheme to solve for evolution of
        a weakness profile eta_i^j and its eroding surface position phi^j.
        """
        W   = self.W
        eta = self.eta_array
        phi = self.phi_array
        chi = self.chi_array
        tau = self.tau_array
        Delta_chi = self.Delta_chi
        Delta_tau = self.Delta_tau
        for j,tau_step in enumerate(tau[:-1]):
            self.j = j
            f = np.int64(phi[j]/Delta_chi)
            f_right = phi[j]/Delta_chi-f
            f_left = 1.0-f_right
            fp1 = f+1
            fp2 = f+2
            if fp2>=self.n_chi_domain:
                break
            self.nu_array[j] = (f_left*eta[j,f]+f_right*eta[j,fp1])
            Delta_phi_j = (self.nu_array[j]*Delta_tau)/(W*2)
            phi[j+1] = phi[j]+Delta_phi_j
            eta[j+1,f:-1] = (
                  eta[j,f:-1] 
                + (Delta_phi_j*(eta[j,fp1:]-eta[j,f:-1]))/(Delta_chi)
                +  Delta_tau*negExpH(chi[f:-1]-phi[j],Delta_chi)
                )
            eta[j+1,-1] = (
                  eta[j,-1] 
                + (Delta_phi_j*(eta[j,-1]-eta[j,-2]))/(Delta_chi)
                +  Delta_tau*negExpH(chi[-1]-phi[j],Delta_chi)
            )
        self.nu_array[j+1] = (f_left*eta[j,f]+f_right*eta[j,fp1])
        di = self.nu_array.shape[0]//10
        self.nu_s_bar = np.mean(self.nu_array[4*di:6*di])
        self.pdict.update({nu_s_bar:self.nu_s_bar})
