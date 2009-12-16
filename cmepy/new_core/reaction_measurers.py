"""
Implementation of measurers for reaction count variables.
"""

import numpy

import cmepy.recorder.measurer as measurer

GROUP_NAME_REACTIONS = 'reactions'

def all_measurers():
    return [ReactionMarginalMeasurer,
            ReactionExpectationMeasurer,
            ReactionStdDevMeasurer,
           ]

class ReactionMarginalMeasurer(measurer.Measurer):    
    def __init__(self, recorder):
        measurer.Measurer.__init__(self,
                                   GROUP_NAME_REACTIONS,
                                   measurer.OUTPUT_NAME_MARGINAL,
                                   recorder)
    
    def compute(self, var_name, t_current, p_current, prereq_measurements):
        # figure out which dimension this variable name corresponds to
        m_info = self.recorder.get_measurement_info(GROUP_NAME_REACTIONS,
                                                    var_name)
        dim = m_info.dim
        # compute marginal distribution for reaction count of dimension dim
        marginal = numpy.array(p_current)
        axis = 0
        for reduce_dim in xrange(self.recorder.dims):
            if reduce_dim == dim:
                axis += 1
                continue
            marginal = numpy.sum(marginal, axis)
        return marginal

class ReactionExpectationMeasurer(measurer.Measurer):    
    def __init__(self, recorder):
        measurer.Measurer.__init__(self,
                                   GROUP_NAME_REACTIONS,
                                   measurer.OUTPUT_NAME_EXPECTATION,
                                   recorder)
    
    def get_prereq_keys(self, var_name):
        marginal_key = (GROUP_NAME_REACTIONS,
                        measurer.OUTPUT_NAME_MARGINAL,
                        var_name)
        return (marginal_key,)
    
    def compute(self, var_name, t_current, p_current, prereq_measurements):
        marginal_key = (GROUP_NAME_REACTIONS,
                        measurer.OUTPUT_NAME_MARGINAL,
                        var_name)
        marginal = prereq_measurements[marginal_key]
        
        # figure out which dimension this variable name corresponds to
        m_info = self.recorder.get_measurement_info(GROUP_NAME_REACTIONS,
                                                    var_name)
        dim = m_info.dim
        
        count_lo = self.recorder.args['norigin'][dim]
        count_hi = count_lo + p_current.shape[dim]
        counts = numpy.arange(count_lo, count_hi)
        
        expectation = numpy.sum(marginal * counts)
        return expectation

class ReactionStdDevMeasurer(measurer.Measurer):
    def __init__(self, recorder):
        measurer.Measurer.__init__(self,
                                   GROUP_NAME_REACTIONS,
                                   measurer.OUTPUT_NAME_STD_DEV,
                                   recorder)
    
    def get_prereq_keys(self, var_name):
        marginal_key = (GROUP_NAME_REACTIONS,
                        measurer.OUTPUT_NAME_MARGINAL,
                        var_name)
        return (marginal_key,)
    
    def compute(self, var_name, t_current, p_current, prereq_measurements):
        marginal_key = (GROUP_NAME_REACTIONS,
                        measurer.OUTPUT_NAME_MARGINAL,
                        var_name)
        marginal = prereq_measurements[marginal_key]
        
        # figure out which dimension this variable name corresponds to
        m_info = self.recorder.get_measurement_info(GROUP_NAME_REACTIONS,
                                                    var_name)
        dim = m_info.dim
        
        count_lo = self.recorder.args['norigin'][dim]
        count_hi = count_lo + p_current.shape[dim]
        counts = numpy.arange(count_lo, count_hi)
        
        expectation = numpy.sum(marginal * counts)
        
        std_dev = numpy.sqrt(numpy.sum(((counts-expectation)**2)*marginal))
        
        return std_dev
