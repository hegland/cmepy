CmePy v0.2
==========
--------------------------------------------
a Chemical Master Equation solver for Python
--------------------------------------------

Features
~~~~~~~~
 *   models can be defined using species or reaction counts
 *   both dense 'rectangular' and sparse state spaces are supported
 *   error due to state space truncation may be tracked with an FSP-style 'sink' state
 *   reaction propensities may be scaled by time dependent coefficients
 *   common statistical results are easily obtained

Dependencies
~~~~~~~~~~~~
CmePy was developed for Python_ 2.5, and depends upon the following packages:

 *   SciPy_ 0.7
 *   Numpy_ 1.2.1
 *   Matplotlib_ (required only to plot example results)

CmePy *should* be compatible with Python_ 2.6, provided SciPy_ 0.7 and
Numpy_ 1.3 are used.

Obtaining CmePy
~~~~~~~~~~~~~~~
Provided you already have the Git_ version control system installed, you may
check out a copy of the latest version of CmePy from GitHub, via::

	git clone git://github.com/fcostin/cmepy.git

Testing and Installation
~~~~~~~~~~~~~~~~~~~~~~~~
Once CmePy has been obtained, the package can be tested and installed
by running the **setup.py** script via Python_ as follows::

    python setup.py test
    python setup.py install

.. _Python: http://www.python.org/
.. _SciPy: http://www.scipy.org/
.. _Numpy: http://numpy.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _Git: http://git-scm.com/