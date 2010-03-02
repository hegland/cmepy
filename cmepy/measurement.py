"""
stores marginal distributions for random variables
"""

class Measurement(object):
    """
    Stores the marginal distributions for a random variable.
    
    Measurement instances are created by cmepy.recorder.CmeRecorder objects.
    
    The random variable is specified by its name and transform when
    instantiating a Measurement object.
    
    Statistics, computed from the marginal distributions, may be accessed
    as attributes (this access is read only and not cached, it is computed
    each time it is accessed).
    """
    def __init__(self, name = None, transform = None):
        """
        Creates a measurement for the specified random variable.
        """
        object.__init__(self)
        self.name = name
        if transform is None:
            self.transform = lambda x : x
        else:
            self.transform = transform
        self.times = []
        self.distributions = []
    
    def __len__(self):
        return len(self.times)
    
    def write(self, t, p):
        """
        Writes the time t and marginal distribution derived from p.
        """
        self.times.append(t)
        self.distributions.append(p.map(self.transform))
    
    def get_statistic(self, stat_name):
        """
        Returns list of statistic values over time.
        
        Alternatively, statistics may be obtained directly as attributes
        of the Measurement instance m, that is,
        
        m.get_statistic(stat_name) <=> m.stat_name
        """
        statistic = []
        for d in self.distributions:
            if stat_name in d.statistics:
                statistic.append(d.statistics[stat_name]())
            else:
                raise KeyError(str(stat_name))
        return statistic
    
    def __getattribute__(self, attrname):
        """
        Allows statistics to be accessed as if they were attributes
        """
        try:
            return object.__getattribute__(self, attrname)
        except AttributeError:
            return self.get_statistic(attrname)
