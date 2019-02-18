"""
---------------------------------------------------------------------

Config to run :py:mod:`wme` in `IPython`_.

Sets up `IPython`_ environment if we're running :py:mod:`wme` 
in a `Jupyter notebook`_ or `Jupyter QtConsole`_. 

 - prepares Matplotlib to display inline and at a 'retina' resolution -- if this
   is not available, a benign error report is made and progress continues
 - enables automatic reloading of :py:mod:`wme` 
   (in case the code has been modded) when 
   a notebook is re-run in-situ



---------------------------------------------------------------------

Requires `matplotlib`_ and `IPython`_.

Uses IPython extensions `autoreload`_.

The  `autoreload`_ extension forces the :py:mod:`wme` 
package to be reloaded on 
restart. This makes code modding and subsequent rerunning of a notebook
smooth and seamless. It is not needed for normal operation, and if unavailable processing 
continues regardless.


---------------------------------------------------------------------

.. _matplotlib: https://matplotlib.org/
.. _autoreload: https://ipython.org/ipython-doc/3/config/extensions/autoreload.html
.. _IPython: https://ipython.readthedocs.io/en/stable/
.. _Jupyter notebook: https://jupyter-notebook.readthedocs.io/en/stable/
.. _Jupyter QtConsole: https://qtconsole.readthedocs.io/en/stable/



"""

# Jupyter `%magic` commands `%load_ext`, `%aimport`, and `%autoreload` 
#  are needed here to force the notebook to reload the `streamline` module, 
#  and its constituent modules, as changes are made to it.
# Force module to reload

# print('Initializing')

import matplotlib as mpl

try:
    get_ipython().magic("config InlineBackend.figure_format = 'retina'")
except NameError as error:
#     print('Error trying to invoke get_ipython(), possibly because not running IPython:', 
#           error)
    pass
except:
    print('Possibly benign error trying to config Matplotlib backend')
    import traceback
    print(traceback.format_exc())
    pass
 
try:
#     get_ipython().magic('matplotlib notebook')
    get_ipython().magic('matplotlib inline')
except NameError as error:
#     print('Error trying to invoke get_ipython(), possibly because not running IPython:', 
#           error)
    pass
except:
    print('Possibly benign error trying to config Matplotlib backend')
    import traceback
    print(traceback.format_exc())
    pass
 
try:
    get_ipython().magic('load_ext autoreload')
    get_ipython().magic('autoreload 2')
    get_ipython().magic('aimport wme')
except NameError as error:
#     print('Error trying to invoke get_ipython(), possibly because not running IPython:', 
#           error)
    pass
except:
    print('Possibly benign error trying to config autoreload')
    import traceback
    print(traceback.format_exc())
    pass

