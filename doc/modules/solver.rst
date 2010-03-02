.. _module-solver:

===================
:mod:`cmepy.solver`
===================
The :mod:`cmepy.solver` module computes numerical solutions to the
Chemical Master Equation (CME).

Creating a solver
~~~~~~~~~~~~~~~~~
Solver instances are obtained by calling the function
:func:`cmepy.solver.create`. This requires the two arguments ``model``
and ``sink``:

 * The ``model`` argument is the model representing the system of reactions
   to solve. To create a model, use :func:`cmepy.model.create`. For more
   information regarding the model format, see :ref:`module-model`.
   
 * The ``sink`` argument is a Boolean flag. If sink is ``True``, the solver
   will include a 'sink' state used to accumulate any probability that may
   flow outside the domain. This can be used to measure the error in the
   solution due to truncation of the domain. If sink is ``False``, the solver
   will not include a 'sink' state, and probability will be artificially
   prevented from flowing outside of the domain.

To specify the initial conditions for the CME, the optional arguments
``p_0`` and ``t_0`` may be passed to :func:`cmepy.solver.create`:

 * ``p_0`` : the initial probability distribution. If given, this must be a
   dictionary from the states in the state space (with respect to the state
   space used by the model) to the initial probabilities. By default, when
   ``p_0`` is not specified, then the argument ``model`` must contain the
   attribute ``model.initial_state``. In this case, the initial probability
   distribution ``p_0`` is defined to be all probability concentrated at the
   initial state, that is, ::
      
      p_0 = { model.initial_state : 1.0 }

 *  ``t_0`` : the initial solution time. By default, ``t_0 = 0.0``.

Finally, two additional optional key word arguments may be passed. These are
``time_dependencies`` and ``domain_states``:

 * ``time_dependencies`` : time dependency data. If given, this must be a
   dictionary of functions taking a time ``t`` and returning a scalar,
   keyed by sets of reaction indices.
   
   See :ref:`time-dependent-propensity-functions`.

 * ``domain_states`` : states to include in the domain. If given, this
   must be an array of states.
   
   See :ref:`sparse-state-spaces`.

Using a solver
~~~~~~~~~~~~~~
Suppose we have obtained a solver instance for the model ``m``, using the
function :func:`cmepy.solver.create`::
    
    from cmepy import model, solver
    
    m = model.create(
        ... # model definition goes here
    )
    
    s = solver.create(m, sink = True)

The solution of the solver ``s`` can be advanced to time ``t`` by calling::

    s.step(t)

Typically, the solver is stepped over a range of times inside a for loop.
For example, to step the solver forward from time 0 to time 10, using 11 evenly
spaced steps (including time 0), use::
    
    import numpy
    
    time_steps = numpy.linspace(0.0, 10.0, 11)
    
    for t in time_steps:
        s.step(t)

The current solution and solution time stored by the solver can be accessed
using the atttributes ``s.y`` and ``s.t`` respectively. If the solver ``s``
was created with the ``sink = True`` flag set, then ``s.y`` will have a value
of the form::

    (p, p_sink) = s.y

where ``p`` is the approximate solution of the CME and ``p_sink`` is the
probability contained by the sink state.
Here, ``p`` is a probability distribution, which is represented as a dictionary
of probabilities, keyed by states in the domain, while ``p_sink`` is a scalar
containing the net probability that has 'leaked' outside the domain into the
sink state.

Conversely, if the solver ``s`` was created with the ``sink`` flag set to
``False``, then ``s.y`` will have a value of the form::

    p = s.y

where ``p`` is a probability distribution of the same form as described above.

Complications
~~~~~~~~~~~~~
For a solver ``s``, the time ``t`` passed to ``s.step(t)`` must not be less than
the solver's current solution time, ``s.t``, otherwise a ``ValueError`` will be
raised.

Additionally, in a rather ugly complication, if the time ``t`` is 'too large',
the underlying 'VODE' ODE solver, obtained from :func:`scipy.integrate.ode`,
may fail to converge. If this occurs, a ``RuntimeError`` will be raised, and
VODE will display an error message. If such errors occur, try reducing the size
of the time steps.
