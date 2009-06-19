"""
Definition of interface to common Solver functionality.
"""


class Solver():
    """
    Interface for common Solver functionality
    """
    def __init__(self, model):
        """
        initialise the Solver with the given model
        """
        assert('propensities' in model), ('model is missing '
            +'a \'propensities\' entry')
        self.model = model
        self.p0 = None
        self.t0 = None
        self.t = None
        # set initial values to defaults...
        self.set_initial_values()
    
    def set_initial_values(self, t0 = 0.0, p0 = None):
        """
        set the initial values
        """
        self.t0 = t0
        self.t = t0
        self.p0 = p0

    def set_solver_params(self, **args):
        """
        set any solver specific parameters
        """
        raise NotImplementedError
    
    def step(self, t):
        """
        advance the solution forward to time t
        """
        raise NotImplementedError
    
    def get_p(self):
        """
        returns a copy of the solution p.
        """
        raise NotImplementedError
    
    def get_t(self):
        """
        returns the current solution time t.
        """
        return self.t
