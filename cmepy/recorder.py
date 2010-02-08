"""
CmeRecorder is a utility class for computation of measurements.
"""

import itertools
from cmepy.statistics import Distribution
from cmepy.measurement import Measurement

def create(*targets):
    """
    alias for CmeRecorder.__init__
    """
    return CmeRecorder(*targets)

class CmeRecorder(dict):
    """
    CmeRecorder is a utility class to compute common measurements,
    such as marginals, expected values, and standard deviations, from
    a given distribution p. 
    """
    
    def __init__(self, *targets):
        """
        Initialise CmeRecorder, optionally using the given sequence of targets.
        
        Each argument should be of the form
            target = (group, variables, transforms)
        where group, variables and transforms are valid arguments
        for the add_target method.
        """
        dict.__init__(self)
        for target in targets:
            if len(target) != 3:
                lament = 'expected (group, variables, transforms), received %s'
                raise ValueError(lament % str(target))
            self.add_target(*target)
        
    def add_target(self, group, variables, transforms=None):
        """
        rec.write(group, variables[, transforms]) : 
        """
        
        if group not in self:
            self[group] = {}
        
        if transforms is None:
            def f(dim):
                return lambda state : state[dim]
            transforms = tuple(f(i) for i in xrange(len(variables)))
        else:
            def pre_star(f):
                return lambda state : f(*state)
            transforms = tuple(pre_star(f) for f in transforms)
        
        if len(variables) != len(transforms):
            raise ValueError('variables and transforms length mismatch')
        
        for var, transform in itertools.izip(variables, transforms):
            self[group][var] = Measurement(var, transform)
        
    def write(self, t, p):
        """
        rec.write(t, p) : records measurements of time t and distribution p
        """
        
        d = Distribution(p)
        
        for group in self:
            for var in self[group]:
                measurement = self[group][var]
                measurement.write(t, d)
    
    def measurements(self, group):
        """
        rec.measurements(group) -> iterator over group measurements
        """
        return self[group].itervalues()

def display_plots(rec, group, statistics = None, title = None):
    """
    plot and display statistics from specified recorder measurement group.
    
    Requires pylab (eg the matplotlib package) to be installed.
    """
    
    if statistics is None:
        statistics = ('expected_value', 'standard_deviation')
    if title is None:
        title = 'Results'
    
    import pylab
    
    for statistic in statistics:
        pylab.figure()
        for measurement in rec.measurements(group):
            pylab.plot(measurement.times,
                       measurement.get_statistic(statistic),
                       label = measurement.name)
        pylab.legend()
        plot_title = ': '.join((title, str(statistic)))
        pylab.title(plot_title)
    pylab.show()
