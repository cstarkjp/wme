Installation
============

 .. _repository: https://github.com/cstarkjp/WeatheringMediatedErosion
 .. _python package: 
       https://github.com/cstarkjp/WeatheringMediatedErosion/tree/master/python
 .. _modules: 
       https://github.com/cstarkjp/WeatheringMediatedErosion/tree/master/python/wm_erosion

First, Git clone the `repository`_:

.. code-block:: none
	
		git clone https://github.com/cstarkjp/WeatheringMediatedErosion.git

Then append the location of 	the `python package`_  to your PYTHONPATH environmental 
variable. For Linux/MacOSX users, mod your .bashrc something like this:

.. code-block:: none
	
		export WMEHOME="${HOME}/<path_to_clone>/WeatheringMediatedErosion"
		export PYTHONPATH="${PYTHONPATH}:${WMEHOME}/python"

That's it.  

If Python fails to import the :mod:`wm_erosion` package, i.e. when doing this:

.. code-block:: python
	
		import wm_erosion as wme
		
it means the PYTHONPATH is not set correctly in the shell from which Python was invoked.
Check the environment variable against the full path to the ``python/`` parent folder of 
the `modules`_.

