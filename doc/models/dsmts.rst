==================================
:mod:`dsmts` : a birth-death model
==================================

Overview
~~~~~~~~
The :mod:`cmepy.models.dsmts` module defines the model
``DSMTS_001_01``. This model is adapted from the
`Discrete Stochastic Models Test Suite
<http://code.google.com/p/dsmts/>`_ .

This 'birth-death' model is defined as the system of reactions:

.. math::

   X & \xrightarrow{0.1} 2 X \; , \\
   X & \xrightarrow{0.11} \star \; . \\

The model is initialised with 100 copies of the species :math:`X`.

Source code
~~~~~~~~~~~
.. literalinclude:: ../../cmepy/models/dsmts.py
