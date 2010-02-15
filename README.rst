CmePy v2.0
==========
a Chemical Master Equation solver
---------------------------------

Features
~~~~~~~~
 * state spaces can be defined using species or reaction counts
 * both dense 'rectangular' and sparse state spaces are supported
 * error due to state space truncation may be tracked with an FSP-style
   'sink' state
 * reaction propensities may be scaled by time dependent coefficients
 * common statistical results are easily obtained (e.g. variance of a
   species count)

Installation
~~~~~~~~~~~~

Dependencies
------------
CmePy was developed with:
* Python 2.5 <http://www.python.org/>
* SciPy 0.7 <http://www.scipy.org/>
* Numpy 1.2.1 <http://numpy.scipy.org/>
* matplotlib 1.2.1 <http://matplotlib.sourceforge.net/>
  (only used by example scripts to plot results)

CmePy *should* be compatible with Python 2.6, provided SciPy 0.7 and Numpy 1.3
are used.

Testing and Installation
------------------------
Once \cmepy{} has been obtained, the package can be tested and installed
by running the **setup.py** script via Python as follows:

    python setup.py test
    python setup.py install
