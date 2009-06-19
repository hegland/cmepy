"""
Implementation of common measurer functionality.
"""

OUTPUT_NAME_EXPECTATION = 'expectation'
OUTPUT_NAME_STD_DEV = 'std_dev'
OUTPUT_NAME_MARGINAL = 'marginal'

class Measurer():
    """
    Base class for the various measurer implementations.
    """
    
    def __init__(self, group_name, output_name, recorder):
        self.recorder = recorder
        self.group_name = group_name
        self.output_name = output_name
    
    def compute(self, var_name, t_current, p_current, preqreq_measurements):
        """
        Outline of interface for compute method of derived classes.
        
        Arguments:
        
        var_name : name of variable to compute
        t_current : the current solution time
        p_current : the current solution CME probability distribution
        prereq_measurements : dictionary of prerequisite measurements,
            keyed by (group_name, output_name, variable_name).
        
        Should return the measurement computed
        """
        raise NotImplementedError
    
    def get_prereq_keys(self, var_name):
        """
        By default, measurer requires no pre-requisite measurements in order
        to compute.
        """
        return []