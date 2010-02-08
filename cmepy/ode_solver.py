import numpy
import scipy.integrate

class Solver(object):
    def __init__(self, dy_dt, y_0, t_0 = 0.0, ode_config_callback = None):
        """
        Initialise a Solver using the supplied derivative function dy_dt,
        initial value y_0, and (optional) initial time t_0.
        
        Optionally, the ode_config_callback argument may also be specified.
        This should be a function of the form
        
            ode_config_callback(ode_instance)
        
        where ode_instance is the scipy.integrate.ode object managed
        internally by the solver.
        """
        self._dy_dt = dy_dt
        self._ode = None
        self._t = t_0
        self._y_0 = y_0
        self._custom_packing = False
        self._ode_config_callback = ode_config_callback
    
    def set_packing(self, pack, unpack, transform_dy_dt=True):
        """
        Convenience routine allowing automatic change of variables.
        
        pack : map from y variables to x variables
        unpack : map from x variables to y variables
        
        x must be a one-dimensional array or scalar, while y is arbitrary
        
        If transform_dy_dt is set to True (the default), the solver will then
        internally solve the IVP of (dx_dt, x_0), defined by
            dx_dt(t, x) := pack(dy_dt(t, unpack(x))
            x_0         := pack(y_0)
        
        Conversely, if transform_dy_dt is set to False, dx_dt will be defined
        directly as dy_dt without packing / unpacking, that is,
            dx_dy(t, x) := dy_dt(t, x)
        
        Accessing the y attribute of the solver will return the current solution
        x transformed via unpack, that is, y = unpack(x).
        """
        assert pack is not None
        assert unpack is not None
        
        if self._ode is not None:
            lament = 'set_packing() must be called prior to step() or accessing solution y'
            raise RuntimeError(lament)
        
        self._custom_packing = True
        self._transform_dy_dt = transform_dy_dt
        self._pack = pack
        self._unpack = unpack
        
        return self
    
    def _initialise_ode(self):
        """
        Internal method, used to initialise scipy.integrate.ode instance when
        necessary.
        """
        
        if self._custom_packing:
            packed_y_0 = self._pack(self._y_0)
            if self._transform_dy_dt:
                def packed_dy_dt(t, y):
                    return self._pack(self._dy_dt(t, self._unpack(y)))
            else:
                packed_dy_dt = self._dy_dt
        else:
            packed_y_0 = self._y_0
            packed_dy_dt = self._dy_dt
            
        ode = scipy.integrate.ode(packed_dy_dt)
        ode.set_integrator('vode', method='bdf')
        ode.set_initial_value(packed_y_0, self._t)
        if self._ode_config_callback is not None:
            self._ode_config_callback(ode)
        self._ode = ode
        self._y_0 = None
    
    @property
    def dy_dt(self):
        """
        Read-only property, returning the differential equations dy_dt.
        """
        return self._dy_dt
    
    @property
    def t(self):
        """
        Read-only property, returning the current solution time t.
        """
        return self._t
    
    @property
    def y(self):
        """
        Read-only property, returning a *copy* of the current solution y.
        """
        if self._ode is None:
            self._initialise_ode()
        # ensure self._y is a *copy* of the solver's current solution
        y = numpy.array(self._ode.y)
        if self._custom_packing:
            y = self._unpack(y)
        return y
    
    def step(self, t):
        """
        Advances the current solution to the time t.
        
        Values of t less that the current solution time are illegal and will
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
        
        self._ode.integrate(t)
        if not self._ode.successful():
            complaint = ('ODE integration failure '+
                         '(look for messages from DVODE / vode)')
            raise RuntimeError, complaint
        self._t = t
