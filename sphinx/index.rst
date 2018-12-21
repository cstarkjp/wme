 .. py:module:: wme

Weathering-mediated erosion
################################################################

Exploration of weathering-driven weakening on rock erosion.


.. image:: images/nu_s.png 
   :target: https://github.com/cstarknyc/wme/blob/master/experiments1d/Results1d

The focus is on environments with no soil accumulation such as on 
steep bedrock channel walls. 
The current treatment is 1d only, and for expository purposes it simplistically assumes 
an exponentially decaying weakening profile with depth into the rock.

The `code`_ is provided as a `Python3 package`_ and an `IPython Jupyter notebook`_.

-----------------------------------------------------------
Set up & execution
-----------------------------------------------------------

 
.. toctree::
   :maxdepth: 1
   
   Prerequisites <documentation/prerequisites>
   Installation <documentation/install>
   Running <documentation/howtorun> 



-----------------------------------------------------------
Code
-----------------------------------------------------------

The :mod:`wme` Python package consists of the following modules:


.. toctree::
   :maxdepth: 1

   modules/initialize   
   modules/symbols
   modules/data
   modules/theory
   modules/solve1d
   modules/solve1p1d
   modules/plot
   modules/save
  
.. only:: html 

* :ref:`genindex` of modules, methods and variables

..  _Python3 package: 
      https://github.com/cstarknyc/wme/tree/master/python
.. _IPython Jupyter notebook: 
      https://github.com/cstarknyc/wme/blob/master/experiments1d/wme1d.ipynb
.. _code: https://github.com/cstarknyc/wme
.. _example results: 
      https://github.com/cstarknyc/wme/blob/master/experiments1d/Results1d


