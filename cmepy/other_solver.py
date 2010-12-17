import numpy
import scipy

#from cmepy import solver
from cmepy.ode_solver import Solver

class SolverOther(Solver):
    """
    SolverOther is a wrapper for several other time solvers, in particular:
    - explicit Euler
    - implicit Euler
    - Heun (explicit two step method)
    - implicit two step method from Deulfhard et. al.
    """
    def __init__(self, dy_dt, y_0, t_0 = 0.0, ode_config_callback = None, **args):
        """
        Initialise a Solver using the supplied derivative function dy_dt,
        initial value y_0, and (optional) initial time t_0.

        !!! NOT REALLY USED !!! UNTESTED !!!
        Optionally, the ode_config_callback argument may also be specified.
        This should be a function of the form

            ode_config_callback(ode_instance)

        where ode_instance is the scipy.integrate.ode object managed
        internally by the solver. If specified, the callback is called
        by the Solver just after ode_instance is instantiated.
        """
        Solver.__init__(self, dy_dt, y_0, t_0, ode_config_callback)
        # later for different SVD-Approximations, maybe move to different setting
        if 'rank' in args:
            self._r = args['rank']
        else:
            self._r = None
        if 'threshold' in args:
            self._threshold = args['threshold']
        else:
            self._threshold = None
        if 'sum' in args:
            self._sum = args['sum']
        else:
            self._sum = None    

        if 'use_reaction_matrix' in args:
            use_reaction_matrix = args['use_reaction_matrix']
        else:
            use_reaction_matrix = False

        if 'use_reaction_matrices' in args:
            use_reaction_matrices = args['use_reaction_matrices']
        else:
            use_reaction_matrices = False

        if 'integrator' in args:
            integrator = args['integrator']
        else:
            integrator = 'imp_euler'

        if (integrator == 'exp_euler'):
            self._do_step = self.exp_euler_step
        elif (integrator == 'imp_euler'):
            self._do_step = self.imp_euler_step
            use_reaction_matrix = True
            self._reaction_matrices = args['reaction_matrices']
        elif (integrator == 'heun'):
            self._do_step = self.heun_step
        elif (integrator == 'imp_two_step'):
            self._do_step = self.imp_two_step
            use_reaction_matrix = True
        else:
            raise ValueError, "Unknown integrator: %s" % integrator

        if use_reaction_matrix:
            # sum the matrices together
            # TODO reaction_matrix = sum(...) ?
            reaction_matrix = None
            for term in args['reaction_matrices']:
                if reaction_matrix is None:
                    reaction_matrix = term
                else:
                    reaction_matrix = reaction_matrix + term
            self._reaction_matrix = reaction_matrix

        if use_reaction_matrices:
            self._reaction_matrices = args['reaction_matrices']

    def exp_euler_step(self, t, y):
        '''
        use explicit Euler
        '''
        h = t - self._t
        y = y + h * self.dy_dt(self._t, y)
        return y, -1


    def heun_step(self, t, y):
        '''
        use heun-method, 
        explicit two-step method of second order and double computational complexity
        y' = f(t,y) is approximated by
               y_n+1 = y_n + h/2 * [ f(t, y_n) + f(t+h, y_n + h*f(t, y_n))].
        '''
        h = t - self._t
        f_t_y = self.dy_dt(self._t, y)
        y_next = y + 0.5*h * ( f_t_y + self.dy_dt(t, y + h*f_t_y) );
        return y_next, -1


    def imp_euler_step(self, t, y):
        '''
        use implicit Euler
        '''
        h = t - self._t
        ident_matrix = scipy.sparse.eye(self._reaction_matrix.shape[0],self._reaction_matrix.shape[1])
        matrix = ident_matrix - h*self._reaction_matrix
        solution = scipy.sparse.linalg.dsolve.spsolve(matrix, y)
        y = solution
        return y, -1


    def imp_two_step(self, t, y, **var):
        '''
        two-step implicit scheme used in Deuflhard et. al., with adaptive time stepping
        '''
        h = t - self._t
        ident_matrix = scipy.sparse.eye(self._reaction_matrix.shape[0],self._reaction_matrix.shape[1])
        matrix = ident_matrix - h*self._reaction_matrix
        rhs = h*self._reaction_matrix*y
        update_0 = scipy.sparse.linalg.dsolve.spsolve(matrix, rhs)
        u_1 = y + update_0
        rhs = -h/2.*self._reaction_matrix*update_0
        update_1 = scipy.sparse.linalg.dsolve.spsolve(matrix, rhs)
        y = u_1 + update_1
        tmp_error = numpy.linalg.norm(update_1)
        if 'TOL' in var:
            TOL = var['TOL']
        else:
            TOL = 0.00001
        if 'sigma' in var:
            sigma = var['sigma']
        else:
            sigma = 0.8
        new_time_step = numpy.sqrt(sigma * TOL/tmp_error)*h
        return y, new_time_step


    def step(self, t, **var):
        """
        Advances the current solution to the time t.

        Values of t less than the current solution time are illegal and will
        raise a ValueError.

        If internal ODE solver errors are detected, a RuntimeError will be
        raised, and additional information may be displayed.
        """
        if self._ode is None:
            self._initialise_ode()

        # guard _ode.integrate(t) against successive calls
        # with non-monotonically increasing values of t
        if t < self._t:
            lament = 'Cannot step backwards to a time t (%f) earlier than current solution time (%f)' % (t, self._t)
            raise ValueError(lament)
        if t == self._t:
            return

        self._ode.y, out = self._do_step(t, self._ode.y, **var)

        self._t = t
        self._y = None
        return out
