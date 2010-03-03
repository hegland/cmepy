===========================
:mod:`munk08` : gene toggle
===========================

Overview
~~~~~~~~
This module provides a model of Gardner's gene toggle system, according to
Munsky and Khammash. This gene toggle model is used as the second example
illustrating CmePy's sparse state space support.

See :ref:`sparse-state-space-gene-toggle-example` for an introduction to this
model.

Running the model
~~~~~~~~~~~~~~~~~
This model is defined by the module :mod:`cmepy.model.munk08`.

To solve the model and display results, simply run the :func:`main` function
inside this model. Explicitly, open the Python interpreter and type:

    >>> from cmepy.models import munk08
    >>> munk08.main()

Source
~~~~~~
.. literalinclude:: ../../cmepy/models/munk08.py
