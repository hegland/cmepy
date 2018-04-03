"""
Creates restorable solvers for the Chemical Master Equation (CME).
"""

import cmepy.solver

def create(model, sink, **solver_args):
    """
    Returns a restorable solver for the CME of the given model.
    
    The arguments are the same as for ``cmepy.solver.create``.
    """
    return RestorableSolver(model, sink, **solver_args)

class RestorableSolver(object):
    """
    CME solver with support for setting and restoring state.
    """
    
    def __init__(self, model, sink, **solver_args):
        """
        Create a RestorableSolver object. 
        """
        self.solver = None
        self.model = model
        self.sink = sink
        self.restore_args = dict(solver_args)
        self.restore()
        self.set_restore_point()
    
    def set_restore_point(self, solver = None):
        """
        Sets a restore point using the current solver state.
        
        Optionally, if the argument solver is given, set a restore
        point using the current state of that solver. This consists
        of the solver's current solution, time, and sink probability,
        if available.
        """
        
        if solver is None:
            solver = self
        
        self.restore_args['t_0'] = solver.t
        if self.sink:
            p, p_sink = solver.y
            self.restore_args['sink_0'] = p_sink
            self.restore_args['p_0'] = p
        else:
            self.restore_args['p_0'] = solver.y
        
    
    def restore(self, **solver_args):
        """
        Restore from a previous restore point fixed via set_restore_point
        
        If set_restore_point was not previously called, restores the
        state of the solver from when it was initialised.
        
        Optionally, one or more keyword arguments may be supplied, which
        then override the arguments used to restore the solver, when
        recreating the solver via cmepy.solver.create .
        """
        restore_args = dict(self.restore_args)
        restore_args.update(solver_args)
        
        self.solver = cmepy.solver.create(
            self.model,
            self.sink,
            **restore_args
        )
    
    def step(self, t):
        """
        Advances the current solution to the time t.
        
        Values of t less that the current solution time are illegal and will
        raise a ValueError.
        
        If internal ODE solver errors are detected, a RuntimeError will be
        raised, and additional information may be displayed.
        """
        self.solver.step(t)
    
    @property
    def restore_point_error(self):
        """
        *Read only* property returning error at restore point.
        
        Raises NotImplementedError if ``sink`` flag is not ``True``.
        """
        if not self.sink:
            raise NotImplementedError('only implemented for sink = True')
        return self.restore_args['sink_0']
    
    @property
    def y(self):
        """
        Read-only property, returning the current solution y.
        """
        return self.solver.y
    
    @property
    def t(self):
        """
        Read-only property, returning the current solution time t.
        """
        return self.solver.t
    
    @property
    def dy_dt(self):
        """
        Read-only property, returning the differential equations dy_dt.
        """
        return self.solver.dy_dt
