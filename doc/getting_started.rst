===============
Getting started
===============

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

Testing and installation
~~~~~~~~~~~~~~~~~~~~~~~~
Once CmePy has been obtained, the package can be tested and installed.

To test CmePy, simply run the **test_all.py** script::

    python test_all.py

If one or more of the tests fail, check to see if all CmePy's
dependencies are correctly installed.

When the **test_all.py** script completes with all tests passed,
then CmePy may be installed by running the **setup.py** script as
follows::

    python setup.py install

You'll probably need to use *sudo* to install CmePy globally::

    sudo python setup.py install

Alternatively, you may wish to use VirtualEnv_ to create an isolated
working environment for CmePy.

.. _Python: http://www.python.org/
.. _SciPy: http://www.scipy.org/
.. _Numpy: http://numpy.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _Git: http://git-scm.com/
.. _VirtualEnv: http://pypi.python.org/pypi/virtualenv
