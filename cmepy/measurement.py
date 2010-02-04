class Measurement(object):
    def __init__(self, name = None, transform = None):
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
        self.times.append(t)
        self.distributions.append(p.map(self.transform))
    
    def get_statistic(self, stat_name):
        """
        m.get_statistic(stat_name) -> list of statistic values over time
        
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
        try:
            return object.__getattribute__(self, attrname)
        except AttributeError:
            return self.get_statistic(attrname)
