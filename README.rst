CmePy v0.3
==========
--------------------------------------------
a Chemical Master Equation solver for Python
--------------------------------------------

Features
~~~~~~~~
 *   models can be defined using species or reaction counts
 *   both dense 'rectangular' and sparse state spaces are supported
 *   error due to state space truncation may be tracked with an FSP-style
     'sink' state
 *   reaction propensities may be scaled by time dependent coefficients
 *   common statistical results are easily obtained

Dependencies
~~~~~~~~~~~~
CmePy was developed for Python_ 2.5, and depends upon the following packages:

 *   SciPy_ 0.7
 *   Numpy_ 1.2.1
 *   Matplotlib_ (required only to plot example results)

CmePy also works with Python_ 2.6, provided SciPy_ 0.7 and Numpy_ 1.3 are used.
The latest version of CmePy also works with Python_ 2.7.

Obtaining CmePy
~~~~~~~~~~~~~~~

If you have the Git_ version control system installed, you can 
check out a copy of CmePy directly from GitHub, via::

	git clone git://github.com/hegland/cmepy.git

Testing and Installation
~~~~~~~~~~~~~~~~~~~~~~~~
Once CmePy has been obtained, the package can be tested by running the
**test_all.py** script via Python_ as follows::

    python test_all.py

CmePy may then be installed via the **setup.py** script::

    python setup.py install

More detailed installation tips are available via the
`online documentation <http://hegland.github.com/cmepy/>`_.


Documentation
~~~~~~~~~~~~~
See http://hegland.github.com/cmepy/

Authors
~~~~~~~

 * Reuben Fletcher-Costin ( reuben dot fletchercostin at gmail dot com )
 * Markus Hegland ( markus dot hegland at anu dot edu dot au )


.. _Python: http://www.python.org/
.. _SciPy: http://www.scipy.org/
.. _Numpy: http://numpy.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _Git: http://git-scm.com/
