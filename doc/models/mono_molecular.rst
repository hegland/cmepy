==============================================
:mod:`mono_molecular` : mono molecular systems
==============================================

Overview
~~~~~~~~

The module :mod:`cmepy.models.mono_molecular` defines models for two
mono-molecular systems of reactions. These systems are useful for testing,
as analytic solutions can be derived. The systems are:

.. math::
   
   A \xrightarrow{} B \xrightarrow{} C

and

.. math::

   A \xrightarrow{} B \; ,
   
   B \xrightarrow{} A \; .

Both systems are initialised with 31 copies of the species
:math:`A` and zero copies of all other species. All reactions have
the rate coefficients equal to 1.

The model for the former system is named ``A2B2C`` while the model
for the latter system is named ``A2B2A``.

For example, the latter model can be imported via::

    from cmepy.models import mono_molecular
    
    model = michaelis_menten.A2B2C


Source
~~~~~~
.. literalinclude:: ../../cmepy/models/mono_molecular.py
