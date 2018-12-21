"""
---------------------------------------------------------------------

Graphing tools to viz behavior of weathering-driven weakening of a 1d rock
mass.

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`numpy`
  -  :mod:`matplotlib.pyplot`
  -  :mod:`matplotlib.ticker`
  -  mpl_toolkits.mplot3d_
  -  :mod:`scipy.optimize`

Imports methods from :mod:`wme`  module :mod:`.solve1d` 

Imports symbols from :mod:`.symbols` module

---------------------------------------------------------------------

.. _mpl_toolkits.mplot3d: https://matplotlib.org/api/toolkits/mplot3d.html#axes3d
.. _`Inoue et al (2017)`: https://doi.org/10.1016/j.geomorph.2017.02.018
.. _`Li et al (2016)`: https://doi.org/10.25103/jestr.093.10

"""


import matplotlib as mpl, matplotlib.pyplot as plt, matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.optimize import curve_fit

from .symbols import *
from .solve1d import eta_chi_tau, nu_s_W
from .data import linear_model, weathering_model

fdict = dict()

def create_figure(fig_name):  
    """
    Initialize a :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
    and add to figures dict.
    
    Create a Matplotlib :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure, 
    append it to the figures dictionary,
    set its size and dpi, set the font size and try to choose the Arial font family.
    
    Args:
        fig_name (str): key for figures dict
    
    Returns:
        :obj:`Matplotlib figure <matplotlib.figure.Figure>`:
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>`  figure
    """
    fig = plt.figure()
    fdict.update({fig_name:fig})
    fig.set_size_inches(7,5)
    fig.set_dpi(100)
    try:
        mpl.rc( 'font', size=12, family='Arial')
    except:
        mpl.rc( 'font', size=12)
    return fig


# def plot_inoue_sigmaT_wetdryN_multiple(fig, ed, text_label=None):
#     """
#     Plot Inoue data on weakness versus proxy time (number of wet/dry cycles).
#     
#     Generate graph of rock weakness versus proxy time inferred from 
#     `Inoue et al (2017)`_
#     data on tensile strength after N wetting and drying cycles.
#     
#     Args:
#         fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
#                   reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>`
#                                         figure 
#         ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
#                                       containing data sets as :mod:`pandas` dataframes
#         text_label (list): text annotation as list of form (x-y coordinate, string, 
#                            font size)
#                            
#     .. _`Inoue et al (2017)`: https://doi.org/10.1016/j.geomorph.2017.02.018
#     
#     """
#     plt.figure(fig.number)
#     
#     df = ed.ddict['inoue']
#     sigmaT  = df.sigmaT
#     wetdryN = df.wetdryN
#     erodibility_sigma2   = df.w_sigma2
#     erodibility_sigma1p5 = df.w_sigma1p5
#     erodibility_sigma2_fit = ed.fdict['inoue'][0]
#         
#     plt.plot(wetdryN,erodibility_sigma2,   label='$n=2$',c='mediumblue', 
#              ls='', marker='o')
#     plt.plot(wetdryN,erodibility_sigma1p5, label='$n=1.5$',
#              color='chocolate',ls='', marker='s')
#     plt.plot(wetdryN,linear_model(wetdryN,*erodibility_sigma2_fit),
#              color='mediumblue',
#              label='$w \sim 1/\\sigma_T^2$')
#     plt.legend(loc='upper left')
#     plt.ylim(0,)
#     plt.xlabel('Proxy time (no. wet/dry cycles)  $N$  [-]')
#     plt.ylabel('Weakness  $w=(\\sigma_T[N]/\\sigma_\mathrm{ref})^{-n}$  [-]')


def plot_inoue_w_wetdryN(fig, ed, text_label=None):
    """
    Plot Inoue et al data on weakness :math:`w` versus 
    proxy time :math:`N` (number of wet/dry cycles).
    
    Generate graph of rock weakness :math:`w` versus proxy time inferred from 
    `Inoue et al (2017)`_
    data on tensile strength :math:`\sigma_T` 
    (normalized by a reference tensile strength :math:`\sigma_\mathrm{ref}`)
    after :math:`N` wetting and drying cycles.
        
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                      containing data sets as :mod:`pandas` dataframes
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)
    
    df = ed.ddict['inoue']
    sigmaT  = df.sigmaT
    wetdryN = df.wetdryN
    erodibility_sigma2   = df.w_sigma2
    erodibility_sigma1p5 = df.w_sigma1p5
    erodibility_sigma2_fit = ed.fdict['inoue'][0]
        
    plt.plot(wetdryN,linear_model(wetdryN,*erodibility_sigma2_fit),
             color='k',
             label='$w \sim 1/\\sigma_T^2$')
    plt.errorbar(np.unique(wetdryN),ed.w_s2_means, label='',
                 xerr=None,
                 yerr=ed.w_s2_stds*1,
                 ecolor='k', mec='w', 
                 color='w', fillstyle='full', 
                 alpha=1,
                 fmt='o', 
                 markersize=0, markeredgewidth=2,
                 elinewidth=1.5,capthick=3,capsize=7)
    plt.plot(np.unique(wetdryN),ed.w_s2_means, label='mean data',
             color='lightgray', fillstyle='full', 
             ls='', marker='o',markeredgecolor='k',ms=15) 
    plt.plot(wetdryN,erodibility_sigma2, label='raw data',  
             color='orange', alpha=0.7,
             ls='', marker='s',markeredgecolor='k',ms=5)

    plt.legend(loc='upper left')
#     plt.ylim(0,)
    plt.xlabel('Proxy time (no. wet/dry cycles)  $N$  [-]')
    plt.ylabel('Weakness  $w=(\\sigma_T[N]/\\sigma_\mathrm{ref})^{-2}$  [-]')

    if text_label is not None:
        axes = plt.gca()
        plt.text(*text_label[0], text_label[1], 
                 color='k', size=text_label[2],
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
        
        
def plot_li_w_wetdryN(fig, ed, text_label=None):
    """
    Plot Li et al  data on weakness :math:`w` versus proxy time :math:`N` 
    (number of wet/dry cycles).
    
    Generate graph of rock weakness :math:`w` versus proxy time :math:`N` 
    inferred from  `Li et al (2016)`_
    data on compressive strength :math:`\sigma_C` 
    (normalized by a reference compressive strength :math:`\sigma_\mathrm{ref}`)
    after :math:`N` wetting and drying cycles
    at a range of confining pressures :math:`P`.
    
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                      containing data sets as :mod:`pandas` dataframes
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)

    df = ed.ddict['li']
    color_list = ('darkblue','darkmagenta','firebrick','red','orange','green')
    marker_list = ('o','s','v','^','8','+')
    n_cols = len(color_list)
    for idx,P__ in enumerate(np.unique(df.P)):
        selection_name = '{0}_{1}_{2}'.format('li','P',P__)
        sigma  = df.sigmaC[df.P==P__]/100
        wetdryN = df.wetdryN[df.P==P__]
        erodibility_sigma   = df.w_sigma2[df.P==P__]
        erodibility_sigma_fit = ed.fdict[selection_name][0]
        
        plt.plot(wetdryN,linear_model(wetdryN,*erodibility_sigma_fit),
                 color=color_list[idx%n_cols],label='')
        plt.errorbar(wetdryN,erodibility_sigma, 
                     label='$P = ${}$\,$MPa'.format(P__),
                     xerr=None,
                     yerr=None,
                     ecolor='k', mec='k', 
                     color=color_list[idx%n_cols], fillstyle='full', 
                     alpha=0.7,
                     fmt=marker_list[idx%n_cols], 
                     markersize=7, markeredgewidth=0.5,
                     elinewidth=1.5,capthick=3,capsize=7)    
    
    
    plt.legend(loc='upper left')
    plt.ylim(0,)
    plt.xlabel('Proxy time (no. wet/dry cycles)  $N$  [-]')
    plt.ylabel('Weakness  $w=(\\sigma_C[N]/\\sigma_\mathrm{ref})^{-2}$  [-]')

    if text_label is not None:
        axes = plt.gca()
        plt.text(*text_label[0], text_label[1], 
                 color='k', size=text_label[2],
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
        
        
def plot_li_w_P(fig, ed, text_label=None):
    """
    Plot Li et al  data on weakness :math:`w` versus proxy depth.
    
    Generate graph of rock weakness :math:`w` versus proxy depth 
    inferred from  `Li et al (2016)`_ data on compressive strength :math:`\sigma_C` 
    at a range of confining pressures :math:`P` after :math:`N` 
    wetting and drying cycles.
    
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                      containing data sets as :mod:`pandas` dataframes
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)
    color_list = ('darkblue','darkmagenta','firebrick','red','orange','green')
    marker_list = ('o','s','v','^','P','p')
    n_cols = len(color_list)
    
    df = ed.ddict['li']
    fit = ed.fdict['li']
    sampled_fit = ed.sdict['li']
    w_ref_vec = np.flipud(sampled_fit[2].T)[:,0]-1
    
    for idx,wetdryN in enumerate(np.flip(np.unique(df.wetdryN))):
        wdN__ = df.wetdryN[df.wetdryN==wetdryN]
        P__ = df.P[df.wetdryN==wetdryN]
        w__ = df.w_sigma2[df.wetdryN==wetdryN]
        P_fit = sampled_fit[1]
        w_fit = np.flipud(sampled_fit[2].T)[idx]
        plt.plot(P_fit, w_fit, color=color_list[idx%n_cols])
        plt.errorbar(P__, w__,
                     label='N = {}'.format(wetdryN),
                     xerr=None,
                     yerr=None,
                     ecolor='k', mec='k', 
                     color=color_list[idx%n_cols], fillstyle='full', 
                     alpha=0.7,
                     fmt=marker_list[idx%n_cols], 
                     markersize=7, markeredgewidth=0.5,
                     elinewidth=1.5,capthick=3,capsize=7)   
        
    plt.legend(loc='upper right')
    plt.ylim(0.,)
    plt.autoscale(enable=True, tight=True, axis='x')
    x_limits = plt.xlim()
    plt.plot(x_limits,(1,1), color='gray', ls=':')
    plt.xlabel('Proxy depth (confining pressure $P$)  [MPa]')
    plt.ylabel('Weakness  $w=(\\sigma_C[N]/\\sigma_\mathrm{ref})^{-2}$  [-]')

    if text_label is not None:
        axes = plt.gca()
        plt.text(*text_label[0], text_label[1], 
                 color='k', size=text_label[2],
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
            
            
def plot_li_w_surface_normed_P(fig, ed, text_label=None):
    """
    Plot Li et al data on weakness :math:`w` (normalized using 2D model) 
    versus proxy depth.
    
    Generate graph of rock weakness versus proxy depth
    inferred from  `Li et al (2016)`_
    data on compressive strength :math:`\sigma_C` 
    at a range of confining pressures :math:`P` after :math:`N` wetting and drying cycles.
    
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                      containing data sets as :mod:`pandas` dataframes
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)
    color_list = ('darkblue','darkmagenta','firebrick','red','orange','green')
    marker_list = ('o','s','v','^','P','p')
    n_cols = len(color_list)
    
    df = ed.ddict['li']
    fit = ed.fdict['li']
    sampled_fit = ed.sdict['li']
    w_ref_vec = np.flipud(sampled_fit[2].T)[:,0]-1
    
    P_fit = sampled_fit[1]
    w_fit = (np.flipud(sampled_fit[2].T)[-1]-1)/w_ref_vec[-1]+1

    plt.errorbar(np.unique(df.P),ed.w_s2normed_means, label='mean data',
                 xerr=None,
                 yerr=ed.w_s2normed_stds*2,
                 ecolor='k', mec='k', 
                 color='lightgray', fillstyle='full', 
                 alpha=0.7,
                 fmt='o', 
                 markersize=14, markeredgewidth=1.5,
                 elinewidth=1.5,capthick=3,capsize=7)
    
    for idx,wetdryN in enumerate(np.flip(np.unique(df.wetdryN))):
        wdN__ = df.wetdryN[df.wetdryN==wetdryN]
        P__ = df.P[df.wetdryN==wetdryN]
        w_normed = df.w_s2normed[df.wetdryN==wetdryN]
        plt.errorbar(P__, w_normed,
                     label='N = {}'.format(wetdryN),
                     xerr=None,
                     yerr=None,
                     ecolor='k', mec='k', 
                     color=color_list[idx%n_cols], fillstyle='full', 
                     alpha=0.7,
                     fmt=marker_list[idx%n_cols], 
                     markersize=7, markeredgewidth=0.5,
                     elinewidth=1.5,capthick=3,capsize=7)   
                
    plt.plot(P_fit, w_fit, color='darkblue', label='$w \sim \exp(-k\chi)$')

    plt.legend(loc='upper right')
    plt.ylim(0.8,)
    plt.autoscale(enable=True, tight=True, axis='x')
    x_limits = plt.xlim()
    plt.plot(x_limits,(1,1), color='gray', ls=':')
    plt.xlabel('Proxy depth (confining pressure $P$)  [MPa]')
    plt.ylabel('Normalized weakness  $w(\\tau,\\chi)/w(0,\\tau)$  [-]')
    
    if text_label is not None:
        axes = plt.gca()
        plt.text(*text_label[0], text_label[1], 
                 color='k', size=text_label[2],
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
        
        
def plot_li_w_wetdryN_P(fig, ed, model_surface, text_label=None):
    """
    Plot Li et al data in 3D on weakness versus proxy depth and confining pressure.
    
    Generate 3D view of 2D surface model & data
    of rock weakness versus proxy depth inferred from 
    `Li et al (2016)`_ data on compressive strength :math:`\sigma_C` 
    at a range of confining pressures :math:`P` after :math:`N` wetting and drying cycles.
    
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        ed (:class:`~.data.ExptData`): instance of experimental :mod:`~.data` class
                                      containing data sets as :mod:`pandas` dataframes
        model_surface (str): 
                key to dict (stored in :attr:`ed`) reference in  to 
                2D  function regressed in 
                :meth:`~.data.ExptData.fit_weathering_model` that models
                joint dependence of weakness on proxy depth and proxy time
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number, projection='3d')
    axes = Axes3D(fig)
    
    df      = ed.ddict['li']
    wdN__   = df.wetdryN
    P__     = df.P
    w__     = df.w_sigma2
    wdN_P__ = np.vstack((wdN__,P__))
    X,Y,Z   = ed.fdict[model_surface]

    axes.scatter(wdN__,P__,w__, color='k', s=40)
    
    colors_list = ('y', 'b')
    colors_array = np.empty(X.shape, dtype=str)
    rf = 3
    for y in range(len(X)):
        for x in range(len(Y)):
            colors_array[x, y] = colors_list[(x//rf+y//rf) % 2]
    axes.plot_surface(X,Y,Z,
                      facecolors=colors_array,  alpha=0.35,
                      linewidth=0, antialiased=True)
    plt.xlim(0,)
    plt.ylim(0,)
    axes.w_xaxis.set_pane_color((1,1,1,1))
    axes.w_yaxis.set_pane_color((1,1,1,1))
    axes.w_zaxis.set_pane_color((1,1,1,1))
    axes.view_init(15, 50)
    axes.view_init(25, 70)
    axes.set_xlabel('Proxy time  $N$  [-]')
    axes.set_ylabel('Proxy depth  $P$ [MPa]')
    axes.set_zlabel('Weakness  $w=$  [-]')


def plot_nu_evolution(fig, ew):
    """
    Plot time-evolution of dimensionless erosion rate :math:`\\nu(\\tau)`.
    
    Generate graph of numerical solutions :math:`j`  of dimensionless rock-surface erosion
    rate :math:`\\nu_i^j` over dimensionless time :math:`\\tau^j`. 
    These solutions are provided in the
    class instance :class:`~.solve1d.ErosionWeathering`.
    Legend-label by weathering number for this instance.
        
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        ew (:class:`~.solve1d.ErosionWeathering`): 
                instance of 1d weathering-mediated erosion model :mod:`~.solve1d` class
    """
    plt.figure(fig.number)
    
    plt.plot(ew.tau_array,ew.nu_array,  color='k', lw=1, label='$W=${}'.format(ew.W))
    plt.legend(loc='center right')
    plt.xlabel('Time  $\\tau$  [-]')
    plt.ylabel('Front speed  $\\partial\\varphi/\\partial\\tau$  [-]')


def plot_eta_evolution(fig, ew, tc=40, nd=2, text_label=None):
    """
    Plot 1d evolution of weathering profile :math:`\eta(\chi,\\tau)` 
    undergoing erosion.
    
    Generate graph of selected solutions :math:`j` 
    of 1d weathering-mediated erosion model
    weakness :math:`\eta_i^j = \eta(\chi_i,\\tau^j)`.
    This series of time slices 
    shows the propagation and development of a steady-state form
    of the weathering depth-profile
    as the rock surface is eroded. Quantities are all dimensionless.
    
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        ew (:class:`~.solve1d.ErosionWeathering`): 
                instance of 1d weathering-mediated erosion model :mod:`~.solve1d` class
        tc (int): 
                :math:`\\tau^j`  slicing 'rate'
        nd (int): 
                number of decimal places in :math:`\\tau`  legend label
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)
    
    chi__ = ew.chi_array
    tau__ = ew.tau_array
    eta__ = ew.eta_array
    j__   = ew.j

    tau_slices1 = np.linspace(0,(ew.tau_n_steps-1)//tc,
                              num=5,endpoint=True,dtype=np.int64)
    tau_slices2 = np.linspace((ew.tau_n_steps-1)//tc*2 ,j__+1,
                              num=5,endpoint=True,dtype=np.int64)
    tau_slices = np.concatenate((tau_slices1,tau_slices2))
    cmap = plt.cm.brg.reversed()
    label='$\\tau$={0:3.'+str(nd)+'f}'
    for idx,tau_slice in enumerate(tau_slices):
        eta_slice = eta__[tau_slice]
        chi_front = chi__[eta_slice==0]
        chi_front = 0 if chi_front.shape[0]==0 else chi_front[-1]
        plt.plot(chi__[chi__>=chi_front],eta_slice[chi__>=chi_front],
                 color=cmap(idx/tau_slices.size),
                 label=label.format(tau__[tau_slice]))
        
    axes = plt.gca()
    plt.xlim(chi__[0],chi__[-1])
    plt.xlim(-(chi__[-1]-chi__[0])/30,chi__[-1]*1.08)
    x_limits = plt.xlim()
    y_limits = plt.ylim()
    plt.ylim(y_limits[0]/3,y_limits[1])
    bbox_props = dict(boxstyle='rarrow,pad=0.3', lw=1.5, 
                      fc='white', ec='k')
    t = axes.text((x_limits[1]-x_limits[0])*0.45, (y_limits[1]-y_limits[0])*0.5, 
                  'erosion', ha='right', va='center', rotation=0, color='k',
                  size=12, bbox=bbox_props)
    
    if text_label is not None:
        plt.text(*text_label[0], text_label[1], 
                 color='k', size=text_label[2],
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)

    plt.legend(loc='upper right')
    plt.xlabel('Distance  $\chi$  [-]')
    plt.ylabel('Weakness  $\eta(\chi,\\tau)$  [-]')


def stability_check(tau,nu):
    """
    Visualize stability of numerical solution.
    
    Check numerical stability of time-stepping by plotting a close-up of
    the erosion front speed over time.
    
    Args:
        tau (numpy.ndarray): time slices :math:`\\tau^j` of numerical solution 
        nu (numpy.ndarray): speed of erosion front :math:`\\nu^j` 
                             at each time slice :math:`\\tau^j` 
    """
    plt.plot(tau,nu,'o-')
    plt.ylim(1.204,1.21);
    plt.xlim(4,4.03);


def plot_etas_steadystate(fig, ew):
    """
    Plot steady-state solution of weakness :math:`\eta_s`.
    
    Graph the numerical solution 
    of the 1d weathering-mediated erosion model
    for weakness :math:`\eta_s(\chi_i | W)` 
    as a function of depth from the rock surface :math:`\chi_i`
    for a given value of the weathering number :math:`W`.
        
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>`
                                        figure 
        ew (:class:`~.solve1d.ErosionWeathering`): 
                instance of 1d weathering-mediated erosion model :mod:`~.solve1d` class
    """
    plt.figure(fig.number)
    
    j__ = (ew.j*3)//4
    i_offset = ew.n_chi_domain//7
    phi__  = ew.phi_array[j__]
    phi0__ = int(phi__/ew.Delta_chi)-i_offset
    chi_s__= ew.chi_array[phi0__:]-ew.chi_array[phi0__+i_offset]
    tau__  = ew.tau_array[-1]

    eta_s_numerical  = ew.eta_array[j__,phi0__:]
    eta_s_analytical = eta_chi_tau(chi_s__,tau__,ew.W)
    
    chi_front = chi_s__[eta_s_numerical==0][-1]
    plt.plot(chi_s__[chi_s__>=chi_front],eta_s_numerical[chi_s__>=chi_front],
             color='k', lw=1, label='numerical')
    plt.plot(chi_s__[chi_s__>=chi_front],eta_s_analytical[chi_s__>=chi_front], 
             color='r', lw=2, label='analytical', ls=(0, (4, 5)) )
    plt.xlim(ew.tau_array[0],ew.tau_array[-1])
    
    axes = plt.gca()
    plt.xlim((chi_s__[0],chi_s__[-1]))
    y_limits = plt.ylim()
    plt.ylim(y_limits[0]/3,y_limits[1])
    bbox_props = dict(boxstyle='rarrow,pad=0.3', lw=1.5, 
                      fc='white', ec='DarkGreen')
    t = axes.text(-0.5, (y_limits[1]-y_limits[0])/2, 
                  'front motion', ha='right', va='center', rotation=0, color='DarkGreen',
                  size=12, bbox=bbox_props)
    
    plt.legend()
    plt.xlabel('Distance relative to front  $\chi_s=\chi-\\varphi_s$  [-]')
    plt.ylabel('Weakness  $\eta_s(\chi_s)$  [-]')
    

def plot_etas_steadystate_set(fig, ew_list, chi_max=8):
    """
    Plot a set of steady-state solutions of weakness :math:`\eta_s`.
    
    Graph a set of numerical solution 
    of the 1d weathering-mediated erosion model
    for weakness :math:`\eta_s(\chi_i | W)` 
    as a function of depth from the rock surface :math:`\chi_i`
    for a set of weathering numbers :math:`W`.
    
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>`
                                        figure 
        ew_list (list): 
                list of instances of 1d weathering-mediated erosion model 
                :mod:`~.solve1d` class :class:`~.solve1d.ErosionWeathering`
    """
    plt.figure(fig.number)
    cmap = plt.cm.brg
    
    chi_min = 0
    for idx,(ew,label) in enumerate(ew_list):
        j__ = (ew.j*3)//4
        i_offset = ew.n_chi_domain//10
        phi__  = ew.phi_array[j__]
        phi0__ = int(phi__/ew.Delta_chi)-i_offset
        chi_s__= ew.chi_array[phi0__:]-ew.chi_array[phi0__+i_offset]
        tau__  = ew.tau_array[-1]
        chi_min = min(chi_min,np.min(chi_s__))
        color=cmap(idx/len(ew_list))
        eta_s_numerical  = ew.eta_array[j__,phi0__:]
        chi_front = chi_s__[eta_s_numerical==0][-1]
        plt.plot(chi_s__[chi_s__>=chi_front],eta_s_numerical[chi_s__>=chi_front],
                 lw=1, color=color, label=label)
        
    axes = plt.gca()
    x_limits = plt.xlim()
    y_limits = plt.ylim()
    plt.ylim(y_limits[0]/3,y_limits[1])
    bbox_props = dict(boxstyle='rarrow,pad=0.3', lw=1.5, 
                      fc='white', ec='DarkGreen')
    t = axes.text(-0.5, (y_limits[1]-y_limits[0])/2, 
                  'front motion', ha='right', va='center', rotation=0, color='DarkGreen',
                  size=12, bbox=bbox_props)
    plt.text(chi_min/3,0.5, 'air', color='k', alpha=0.7, size=14, 
             verticalalignment='center', horizontalalignment='right')
    plt.text(-chi_min/3,0.5, 'rock', color='k', alpha=0.7, size=14,
             verticalalignment='center', horizontalalignment='left')
        
    plt.xlim((chi_min,chi_max))
    plt.legend()
    plt.xlabel('Distance relative to front  $\chi_s=\chi-\\varphi_s$  [-]')
    plt.ylabel('Weakness  $\eta_s(\chi_s)$  [-]')


def plot_nus_W(fig, em, do_loglog=True, nus_solns_list=None, text_label=None):
    """
    Plot the 1d model steady-state erosion rate :math:`\\nu_s` 
    versus weathering number :math:`W`.
    
    Graph the functional dependence of dimensionless steady-state 
    erosion rate :math:`\\nu_s` as a function of versus weathering number :math:`W`
    for the 1d weathering-mediated erosion model.
    The analytical solution is plotted as a black curve; numerical solutions are
    plotted as black circles; asymptotic behavior for low and high :math:`W`
    are shown as dashed lines. Explanatory annotations are included.
    
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        em (:class:`~.theory.WeatheringMediatedErosion`): 
                instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        do_loglog (bool): 
            flag whether to use log scales on both axes
        nus_solns_list (list): 
            set of numerical solutions of :math:`\\nu_s`
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)
     
    n_W_pts = 200
    W_array   = np.exp(np.linspace(np.log(0.01),np.log(50),n_W_pts))
    nus_array = np.array([em.nus_eqn_W.rhs.subs({W:W__}) for W__ in W_array])
    y_limits = (nus_array[0]*0.95,nus_array[-1])
 
    plt.plot(W_array, nus_array, color='k', lw=1.5, label='analytical')
    plt.autoscale(enable=True, tight=True)
    axes = plt.gca()
    if do_loglog:
        axes.set_xscale("log", nonposx='clip')
        axes.set_yscale("log", nonposy='clip')
        axes.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
        axes.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%0.0f'))
        plt.plot((0.25,0.25),(y_limits[0],y_limits[1]/3), color='gray', ls=':')
        plt.plot((2.63,2.63),(y_limits[0],y_limits[1]), color='gray', ls=':')
#         plt.ylim(*y_limits)
        plt.text(0.2,0.35, 'low $W$', color='brown',
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
        plt.text(0.2,0.25, '$\\nu_s \\approx 1 + W$', color='brown',
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
        plt.text(0.2,0.15, '$v_s \\approx v_0 + \\dfrac{w_0}{k}$', color='brown',
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
        plt.text(0.52,0.43, 'transitional $W$', color='gray',
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
        plt.text(0.75,0.7, 'high $W$', color='blue',
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
        plt.text(0.78,0.8, '$\\nu_s \\approx \\dfrac{1}{2}+\sqrt{W}$', 
                 color='blue',
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
        plt.text(0.81,0.92, '$v_s \\approx \\dfrac{v_0}{2}+\sqrt{\\dfrac{v_0 w_0}{k}}$', 
                 color='blue',
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
         
    if nus_solns_list is not None:
        for idx,nus_soln in enumerate(nus_solns_list):
            plt.plot(nus_soln.W,nus_soln.nu_s,'o',c='k',
                    label=('numerical' if idx==0 else None))
     
    plt.plot(W_array[W_array<0.7],1+W_array[W_array<0.7], 
             label='low W approx',  ls='--',c='brown')
    plt.plot(W_array[W_array>1],0.5+np.sqrt(W_array[W_array>1]), 
             label='high W approx', ls='--',c='blue')
    if text_label is not None:
        plt.text(0.82,0.25, text_label, 
                 color='k', size=14,
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
     
    plt.legend(loc='upper left')
    plt.xlabel('Weathering number  $W$  [-]')
    plt.ylabel('Erosion rate  $\\nu_s$  [-]')


def plot_nus_W_transition(fig, em, text_label=None):
    """
    Plot the steady-state erosion rate :math:`\\nu_s` relative to its asymptotic behavior.
    
    Visualize the transitional behavior of the steady-state erosion rate :math:`\\nu_s`
    at intermediate weathering numbers :math:`W` by plotting how asymptotes at
    low and high :math:`W` deviate from the full model.
    
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        em (:class:`~.theory.WeatheringMediatedErosion`): 
                instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)
    
    n_W_pts = 200
    W_array    = np.exp(np.linspace(np.log(0.015),np.log(50),n_W_pts))
    nus_array = np.array([em.nus_eqn_W.rhs.subs({W:W__}) for W__ in W_array])
    plt.plot(W_array, 1+0*W_array, color='k', lw=1.5)
    axes = plt.gca()
    plt.plot(W_array[W_array<0.48], 
         nus_array[W_array<0.48]/(1+(W_array[W_array<0.48])),
             color='brown', ls='-', lw=1.5, label='low $W$ approx')
    plt.plot(W_array[W_array>1], 
             nus_array[W_array>1]/(0.5+np.sqrt(W_array[W_array>1])),
             color='b', ls='-', lw=1.5, label='high $W$ approx')
    plt.text(0.18,0.18, 'low $W$', color='brown',
             verticalalignment='center', horizontalalignment='center',
             transform=axes.transAxes)
    plt.text(0.18,0.26, '$\\nu_s \\approx 1 + W$', color='brown',
             verticalalignment='center', horizontalalignment='center',
             transform=axes.transAxes)
    plt.text(0.5,0.53, 'transitional $W$', color='gray',
             verticalalignment='center', horizontalalignment='center',
             transform=axes.transAxes)
    plt.text(0.82,0.85, 'high $W$', color='blue',
             verticalalignment='center', horizontalalignment='center',
             transform=axes.transAxes)
    plt.text(0.82,0.75, '$\\nu_s \\approx \\dfrac{1}{2}+\sqrt{W}$', 
             color='blue',
             verticalalignment='center', horizontalalignment='center',
             transform=axes.transAxes)
    if text_label is not None:
        plt.text(0.82,0.25, text_label, 
                 color='k', size=14,
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
    plt.autoscale(enable=True, tight=True)
    axes.set_xscale('log', nonposx='clip')
    y_limits = (0.9,1.1)
    plt.ylim(*y_limits)
    axes.set_yticks(np.linspace(*y_limits,5))
    axes.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    plt.plot((0.25,0.25),y_limits, color='gray', ls=':')
    plt.plot((2.63,2.63),y_limits, color='gray', ls=':')
    plt.legend(loc='upper left')
    plt.xlabel('Weathering number  $W$  [-]')
    plt.ylabel('Approx erosion rate deviation  $\\nu_s^\mathrm{apx}/\\nu_s$  [-]')


def plot_v0_vs_w0(fig, em, k__=1, text_label=None):
    """
    Plot baseline erosion rate versus baseline weathering rate.
        
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        em (:class:`~.theory.WeatheringMediatedErosion`): 
                instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        k__ (float): 
            weathering-weakening rate profile reciprocal e-folding depth
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)
    
    w0_array = 10**np.linspace(-3,+1,100)
    for vs__ in [0.1,1,3]:
        v0_array = np.array(
            [sy.N(em.v0_eqn_vs_w0.rhs.subs({v_s:vs__, w_0:w0__,k:k__})) 
                                for w0__ in w0_array])
        plt.plot(v0_array,w0_array, label='$v_s=${}'.format(vs__))
    plt.ylabel('$w_0$')
    plt.xlabel('$v_0$')
    plt.legend(loc='lower right')


def plot_v0_vs_etas0(fig, em, text_label=None):
    """
    Plot baseline erosion rate versus surface weakness at steady-state.
        
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        em (:class:`~.theory.WeatheringMediatedErosion`): 
                instance of 1d weathering-mediated erosion theory :mod:`~.theory` class
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)
    
    etas0_array = 10**np.linspace(-2,+1,100)
    for vs__ in [0.5,1,1.5]:
        v0_array = np.array(
            [sy.N(em.v0_eqn_etas0_vs.rhs.subs({v_s:vs__,eta_s0:etas0__})) 
                                for etas0__ in etas0_array])
        plt.plot(etas0_array,v0_array, label='$v_s=${}'.format(vs__))
    plt.xlabel('Surface weakness (degree of weathering)  $\\eta_{s0}$')
    plt.ylabel('Baseline (potential) erosion rate  $v_0$')
    plt.gca().invert_yaxis()
    plt.xlim(0,3)
    plt.ylim(5,0)
    plt.legend(loc='lower right')


def plot_channel_generic(fig, zy_list, text_labels=None, do_equal_aspect=False):
    """
    Plot numerical solutions applied to channel cross-section model (vertical profiles).
        
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure
        zy_list (list): 
            set of numerical solutions to plot
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
        do_equal_aspect (bool): 
            flag whether to use force equal sizing of x and y axis scales
    """
    plt.figure(fig.number)
    
    zy=zy_list[0]
    plt.plot(zy[2], zy[0], label=zy[4], color='k')
    plt.ylabel(zy[1])
    plt.xlabel(zy[3])
    axes = plt.gca()
    if do_equal_aspect:
        axes.set_aspect('equal')
    else:
        pass
    if text_labels is not None:
        for text_label in text_labels:
            plt.text(*text_label[0], text_label[1], 
                     color=text_label[3], size=text_label[2],
                     verticalalignment='center', horizontalalignment='center',
                     transform=axes.transAxes, rotation=text_label[4])
    plt.grid('on',ls=':')
    if len(zy_list)>=2:
        zy=zy_list[1]
        plt.plot(0,0, label=zy[4], color='forestgreen')
    plt.legend()
    
    if len(zy_list)>=2:
        zy=zy_list[1]
        alt_axes = axes.twiny()
        alt_axes.plot(zy[2], zy[0], label=zy[4], color='forestgreen')
        alt_axes.set_xlabel(zy[3], color='forestgreen')


def plot_channel_w0_v0_W(fig, cw, text_label=None):
    """
    Plot numerical solutions applied to channel cross-section model (vertical profiles).
        
    Args:
        fig (:obj:`Matplotlib figure <matplotlib.figure.Figure>`): 
            reference to :mod:`MatPlotLib/Pyplot <matplotlib.pyplot>` figure 
        cw (:class:`~.solve1p1d.ChannelWall`): instance of :mod:`~.solve1p1d` model 
                             class that simulates channel cross-sectional geometry
        text_label (list): text annotation as list of form (x-y coordinate, string, 
                           font size)
    """
    plt.figure(fig.number)

    plt.plot(cw.w0_array/cw.pdict[k], cw.z_array, label='$w_0/k$')
    plt.plot(cw.v0_array, cw.z_array, label='$v_0$')
    x_limits = plt.xlim()
    y_limits = plt.ylim()
    plt.plot([-1,-1], label='$W=w_0/v_0 k$', color='forestgreen')  # dummy
    plt.plot(cw.vs_array, cw.z_array, label='$v_s$', color='k', lw=2)
    plt.xlim(x_limits)
    plt.ylim(y_limits)
    plt.xlabel('Speeds $w_0(z)/k$, $v_0(z)$, $v_s(z)$')
    plt.ylabel('Height above bed  $z$')
    plt.legend(loc='upper center')
    plt.grid('on',ls=':')

    axes = plt.gca()
    alt_axes = axes.twiny()
    alt_axes.plot(cw.W_array,  cw.z_array, label='$W$', color='forestgreen')
    alt_axes.set_xlabel('Weathering number  $W(z)$', color='forestgreen')
    x_limits = axes.get_xlim()
    axes.set_xlim(x_limits[0],x_limits[1]*1.05)

    if text_label is not None:
        plt.text(*text_label[0], text_label[1], 
                 color='k', size=text_label[2],
                 verticalalignment='center', horizontalalignment='center',
                 transform=axes.transAxes)
