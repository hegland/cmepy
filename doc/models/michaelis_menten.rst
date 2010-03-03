===============================================
:mod:`michaelis_menten` : an enzymatic reaction
===============================================

Overview
~~~~~~~~
The :mod:`cmepy.models.michaelis_menten` module defines a model for the 
following simple enzymatic reaction system:

.. math::
   
   E + S \xrightarrow{k_1} C \; ,
   C \xrightarrow{k_2} E + S \; ,
   C \xrightarrow{k_3} E + D \; .

The rate coefficients are defined to be
:math:`k_1 = 0.01`,
:math:`k_1 = 35.0`
and
:math:`k_1 = 30.0`,
while the default initial counts are 50 copies of the species :math:`S` and
10 copies of the species :math:`E`, with all other initial species counts set to
zero.

This model can be used in CmePy as follows::

    from cmepy.models import michaelis_menten
    
    model = michaelis_menten.create_model_michaelis_menten()

This is roughly the same model used in the enzyme kinetics example.
See :ref:`example-enzyme-kinetics` for a detailed example explaining how this
similar model can be defined and solved.

Source
~~~~~~
.. literalinclude:: ../../cmepy/models/michaelis_menten.py
