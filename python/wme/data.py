"""
---------------------------------------------------------------------

Module to import and model Inoue & Li experimental data.

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`numpy`
  -  :mod:`os`
  -  :mod:`pandas`
  -  :mod:`scipy.optimize`

---------------------------------------------------------------------

.. _`Inoue et al (2017)`: https://doi.org/10.1016/j.geomorph.2017.02.018
.. _`Li et al (2016)`: https://doi.org/10.25103/jestr.093.10

"""

import os, pandas as pd, numpy as np
from scipy.optimize import curve_fit

def linear_model(x,m,c):
    """
    Simple linear model of form: :math:`y = m x + c`
    
    Args:
        x (float or numpy.ndarray) : coordinate
        m (float) : gradient
        c (float) : intercept

    Returns:
        float or numpy.ndarray: y
    """
    return m*x+c

def exponential_decay_model(x,m,c):
    """
    Shifted exponential decay model of form: :math:`y = 1 + c \exp(-x/m)`
    
    Args:
        x (float or numpy.ndarray) : coordinate
        m (float) : e-folding scale
        c (float) : magnitude

    Returns:
        float or numpy.ndarray: y
    """
    return 1 + c*np.exp(-x/m)

def weathering_model(wetdryN_P,k,w0,tau0):
    """
    Shifted exponential decay weathering model: 
        :math:`w = 1 + w_0(\\tau+\\tau_0)\exp(-k\chi)`
    
    Args:
        wetdryN_P (numpy.ndarray) : pair :math:`(\\tau,\chi)`
        k (float) : reciprocal e-folding scale :math:`k`
        w0 (float) : magnitude :math:`w_0`
        tau0 (float) : time offset :math:`\\tau_0`

    Returns:
        float: y
    """
    tau = wetdryN_P[0]
    chi = wetdryN_P[1]
    return 1 + w0*(tau+tau0)*np.exp(-k*chi)

class ExptData:
    """
    Rock weathering experimental data import and modeling
    
    Attributes:
        ddict (dict) : data dictionary, to contain experimental results of 
            `Inoue et al (2017)`_ and `Li et al (2016)`_ variously on 
            rock tensile strength, rock compressive strength,
            erodibility, number of wetting/drying cycles, and confining pressure
        fdict (dict) : 1d model dictionary, to contain regression fits of
            weathering model(s) to experimental data
        sdict (dict) : 2d model dictionary, to contain regression fits of
            2d weathering model to experimental data

    """
    def __init__(self):
        """
        Initialize class instance.
        """
        self.ddict = dict()
        self.fdict = dict()
        self.sdict = dict()
        
    def read_excel(self, data_set,
                   dir_name=('..','data'), file_name='Inoue_wetdryN_sigmaT',
                   header=0, skiprows=[1]):
        """
        TBD
        
        Attributes:
            xx (xx) : xxxx
    
        """
        dir_name = os.path.join(*dir_name)
        if not os.path.exists(dir_name):
            print('Cannot find data directory')
            raise
        try:
            df = pd.read_excel(os.path.join(dir_name,file_name+'.xlsx'),
                               header=header, skiprows=skiprows)
        except OSError:  
            print('Cannot find data directory')
            raise
        except:  
            raise
        self.ddict.update({data_set:df})

    def fit_linear_model(self, data_set, x_name, y_name, select=None):
        """
        Regress a 1d model against experimental data.
        
        Perform a linear regression fit of a linear model to the given
        experimental data, such as modeling the degree of rock weakness
        as a linear function of the number of wetting and drying cycles.
        
        Args:
            data_set (str) : xxx
            x_name (numpy.ndarray) : xxx
            y_name (numpy.ndarray) : xxx
            select (int) : xxx
    
        Attributes:
            fdict[data_set] or fdict[selection_name]  (dict element) : xxxx
            w_s2_means (xx) : xxxx
            w_s2_stds (xx) : xxxx
    
        """
        df = self.ddict[data_set]
        if select is not None:
            for selection in np.unique(self.ddict[data_set][select]):
                selection_name = '{0}_{1}_{2}'.format(data_set,select,selection)
                x = df.loc[df[select]==selection][x_name]
                y = df.loc[df[select]==selection][y_name]
                self.fdict[selection_name] = curve_fit(linear_model, x, y)
        else:
            x = df[x_name]
            y = df[y_name]
            self.fdict[data_set] = curve_fit(linear_model, x, y)
            
        self.w_s2_means = self.ddict[data_set].groupby('wetdryN').mean()['w_sigma2']
        self.w_s2_stds  = self.ddict[data_set].groupby('wetdryN').std()[ 'w_sigma2']

            
    def fit_weathering_model(self, data_set, z):
        """
        Regress a 2d model against experimental data.

        
        Args:
            data_set (str) : xxx
            z (numpy.ndarray) : xxx

        Attributes:
            fdict[data_set] (dict element) : xxxx
            sdict[data_set] (dict element) : xxxx
            w_s2normed_means (xx) : xxxx
            w_s2normed_stds (xx) : xxxx
    
        """
        df = self.ddict[data_set]
        wdN_vec  = df.wetdryN
        P_vec    = df.P
        sig_vec  = df.sigmaC/180
        w_vec    = sig_vec**(-2)
        wdN_P_array = np.vstack((wdN_vec,P_vec))
        model_fit = curve_fit(weathering_model, wdN_P_array, w_vec)
        self.fdict[data_set] = model_fit
    
        n_pts = 30
        X = np.linspace(0, wdN_vec.max()*1.1, n_pts)
        Y = np.linspace(0, P_vec.max()*1.1,   n_pts)
        X,Y = np.meshgrid(X, Y)
        X_Y_array = np.vstack((X.reshape(n_pts**2), Y.reshape(n_pts**2)))
        Z = weathering_model(X_Y_array,*model_fit[0]).reshape(n_pts,n_pts)
        self.fdict[data_set+'_'+z+'_surface'] = (X,Y,Z)
        
        P_vec = np.linspace(0,P_vec.max()*1.3)
        X,Y = np.meshgrid(np.unique(wdN_vec), P_vec)
        X_Y_array = np.vstack((X.reshape(X.shape[0]*Y.shape[1]), 
                               Y.reshape(X.shape[0]*Y.shape[1])))
        self.sdict[data_set] \
            = X[0], Y[:,0], \
                weathering_model(X_Y_array,*model_fit[0]).reshape(X.shape[0],Y.shape[1])
        
        pd.options.mode.chained_assignment = None
        w_ref_vec = (self.sdict[data_set][2].T.copy())[:,0]-1
#         w_ref_vec = np.flipud(self.sdict[data_set][2].T)[:,0]-1
        df['w_s2normed'] = 0
        for idx,wetdryN in enumerate(np.unique(df.wetdryN)):
            w__ = df.w_sigma2[df.wetdryN==wetdryN].copy()
            w_normed = (w__-1)/w_ref_vec[idx]+1
            df.w_s2normed[(idx*3):(idx*3+3)] = w_normed

        self.w_s2normed_means = self.ddict['li'].groupby('P').mean()['w_s2normed']
        self.w_s2normed_stds  = self.ddict['li'].groupby('P').std()['w_s2normed']



