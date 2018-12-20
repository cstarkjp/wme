"""
---------------------------------------------------------------------

Mathematical symbols for parameters. 
Make for pretty printing of parameter dict contents in Jupyter.

---------------------------------------------------------------------

Requires Python packages/modules:
  -  :mod:`sympy`

---------------------------------------------------------------------

"""

import sympy as sy

W, w_0, v_0, k, nu_s, v_s = sy.symbols('W w_0 v_0 k \\nu_s v_s', positive=True)
eta, eta_s, eta_s0 = sy.symbols('\\eta \\eta_s \\eta_{s0}', positive=True)
v_r, h, z, w_r, v_b = sy.symbols('v_r h z w_r v_b', positive=True)
v_smooth, v_off = sy.symbols('v_\mathrm{smooth} v_\mathrm{off}', positive=True)
w_smooth, w_off = sy.symbols('w_\mathrm{smooth} w_\mathrm{off}', positive=True)
nu_s_bar = sy.symbols('\\overline{\\nu}_s', positive=True)
Delta_chi = sy.symbols('\\Delta\\chi', positive=True)
Delta_tau = sy.symbols('\\Delta\\tau', positive=True)
chi_domain_size = sy.symbols('\chi_\mathrm{domain}', positive=True)
tau_domain_size = sy.symbols('\\tau_\mathrm{domain}', positive=True)
