How to run
###########

The :mod:`wm_erosion` Python package can be invoked in several ways.


Interactive IPython/Jupyter notebook
------------------------------------------------------------------------

The recommended approach is to deploy :mod:`wm_erosion`  in a Jupyter (browser) 
session and using a Jupyter/IPython notebook. 

Solution of the 1d problem is carried out in the `WeatheringMediatedErosion1d.ipynb`_
available on the `Github repo experiments section`_.


  
In a IPython/Jupyter QtConsole (inline graphics): interactive or non-interactive
--------------------------------------------------------------------------------

Computation of a notebook using :mod:`wm_erosion` can be invoked in 
a Jupyter `QtConsole`_ running an `IPython`_ kernel. 

::

	Jupyter QtConsole 4.3.1
	Python 3.6.4 (default, Dec 21 2017, 20:33:17) 
	Type 'copyright', 'credits' or 'license' for more information
	IPython 6.2.1 -- An enhanced Interactive Python. Type '?' for help.
	
	run WeatheringMediatedErosion1d.ipynb
	
	[... graphs plotted after a few seconds of computation]

Graphical output will (depending on :mod:`initialize <wm_erosion.initialize>` 
settings) be displayed inline.



In a IPython/Jupyter console (external viewer)  
----------------------------------------------------------------------

Similarly, computation of a notebook using :mod:`wm_erosion` can be invoked from 
a Jupyter console running IPython. 

::

	% jupyter-console-3.6 WeatheringMediatedErosion1d.ipynb 
	Jupyter console 5.2.0
	
	Python 3.6.4 (default, Dec 21 2017, 20:33:17) 
	Type 'copyright', 'credits' or 'license' for more information
	IPython 6.2.1 -- An enhanced Interactive Python. Type '?' for help.
	
	In [1]: run WeatheringMediatedErosion1d.ipynb
	
	[... graphs plotted after a few seconds of computation]
	

Graphical output will be pushed to a viewer external to the shell.





.. _Github repo experiments section: 
      https://github.com/cstarknyc/WeatheringMediatedErosion/tree/master/experiments1d
.. _WeatheringMediatedErosion1d.ipynb: 
      https://github.com/cstarknyc/WeatheringMediatedErosion/blob/master/experiments1d/WeatheringMediatedErosion1d.ipynb
.. _WeatheringMediatedErosion: https://github.com/cstarknyc/WeatheringMediatedErosion
.. _QtConsole: https://ipython.org/ipython-doc/3/interactive/qtconsole.html
.. _IPython: http://ipython.org/ipython-doc/3/interactive/

