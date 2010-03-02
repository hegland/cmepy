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
Testing CmePy
-------------
To test CmePy, simply run the **test_all.py** script::

    python test_all.py

If one or more of the tests fail, check to see if all CmePy's
dependencies are correctly installed.

Installing CmePy globally
-------------------------
When the **test_all.py** script completes with all tests passed,
then CmePy may be installed by running the **setup.py** script as
follows::

    python setup.py install

On linux, you will most likely need to use *sudo* to install CmePy globally::

    sudo python setup.py install

If it is not possible or not desirable to install CmePy as a global package,
there are a couple of ways to install CmePy locally. The first way is perhaps
quicker to get started, but is a bit ugly, as it relies on messing about
with environment variables,. The second way is to use the VirtualEnv
package to create and manage a local, isolated working environment for CmePy.

Local installation via environment variables (linux)
----------------------------------------------------
If it is not possible to install CmePy to the Python's global
``site_packages`` directory, CmePy may instead be installed locally to a user's
home directory. The following instructions assume you are using the bash shell.

First, ensure that the path ``${HOME}/lib/python`` exists,
and that the ``${PYTHONPATH}`` environment variable includes
``${HOME}/lib/python``. Then, CmePy may be installed locally as follows::

    python setup.py install --home=${HOME}

Local installation using the VirtualEnv package
-----------------------------------------------
The VirtualEnv package is a tool to create isolated Python environments.
This addresses the far more general problem of managing Python environments
for multiple projects. VirtualEnv can be combined with Pip, Distribute,
and VirtualEnvWrapper to create a very nice system for managing
Python packages.

To learn more about these useful packages, visit their pages on the Python
package index:

 * `Distribute <http://pypi.python.org/pypi/distribute>`_
 * `Pip <http://pypi.python.org/pypi/pip>`_
 * `VirtualEnv <http://pypi.python.org/pypi/virtualenv>`_
 * `VirtualEnvWrapper <http://pypi.python.org/pypi/virtualenvwrapper>`_

.. _Python: http://www.python.org/
.. _SciPy: http://www.scipy.org/
.. _Numpy: http://numpy.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _Git: http://git-scm.com/
