.. _fsp-api-documentation:

==========
:mod:`fsp`
==========

Overview
--------
The sub-package :mod:`cmepy.fsp` contains CmePy's
implementation of the Finite State Projection (FSP) algorithm.
This sub-package contains the following modules:

.. toctree::
   :glob:
      
   fsp/*

Example scripts
---------------
Some scripts containing examples of using CmePy's FSP
implementation to solve the burr08 model are included in
the examples directory.

The example scripts are named:

* examples/fsp_example_1_simple_expansion.py
* examples/fsp_example_2_support_exansion.py

These scripts all use a common function to plot output,
which is defined by the file:

 * fsp_example_util.py

To run the scripts, simply execute them using the Python interpreter. For example::

    python fsp_example_1_simple_expansion.py

Each script will display two plots via matplotlib when complete. The first
plot shows the solution of the CME at a number of times. The second plot
shows the domain used by the FSP solver for these times.

.. Note::
   
   If these example scripts are unable to locate CmePy,
   check that CmePy is correctly installed.

Example script documentation
----------------------------

.. toctree::
   :glob:
      
   ../fsp_examples/*
