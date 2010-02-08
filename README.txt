CmePy release v 2.0
-------------------

CmePy solves the chemical master equation, using SciPy's 'VODE' ODE solver.

Features
--------
 * state spaces can be defined using species or reaction counts
 * both dense 'rectangular' and sparse state spaces are supported
 * error due to state space truncation may be tracked with an FSP-style 'sink' state
 * reaction propensities may be scaled by time dependent coefficients
 * common statistical results are easily obtained (e.g. variance of a species count)

Installation & Usage
--------------------

Please refer to doc/tutorial.pdf for a guide to installation, and an introduction to usage.
Example CmePy scripts are provided inside the examples directory.

