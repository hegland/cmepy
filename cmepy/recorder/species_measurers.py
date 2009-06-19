"""
Implementation of measurers for species count variables.
"""

import numpy
import itertools
import cmepy.recorder.measurer as measurer

GROUP_NAME_SPECIES = 'species'
OUTPUT_NAME_SPARSE_MARGINAL = 'sparse_marginal'

def all_measurers():
    return [SpeciesSparseMarginalMeasurer,
            SpeciesExpectationMeasurer,
            SpeciesStdDevMeasurer,
           ]

class SpeciesSparseMarginalMeasurer(measurer.Measurer):    
    def __init__(self, recorder):
        measurer.Measurer.__init__(self,
                                   GROUP_NAME_SPECIES,
                                   OUTPUT_NAME_SPARSE_MARGINAL,
                                   recorder)
    
    def compute(self, var_name, t_current, p_current, prereq_measurements):
        # figure out which dimension this variable name corresponds to
        m_info = self.recorder.get_measurement_info(GROUP_NAME_SPECIES,
                                                    var_name)
        dim = m_info.dim
        
        indices = numpy.indices(p_current.shape)
        
        species_count_funcs = self.recorder.solver.model['species counts']
        species_count_func = species_count_funcs[dim]
        
        marginal = {}
        indices_flat = [numpy.ravel(ind) for ind in indices]
        for reaction_state in itertools.izip(*indices_flat):
            prob = p_current[reaction_state]
            if prob != 0.0:
                count = species_count_func(*reaction_state)    
                marginal[count] = marginal.get(count, 0.0) + prob
        return marginal

class SpeciesExpectationMeasurer(measurer.Measurer): 
    def __init__(self, recorder):
        measurer.Measurer.__init__(self,
                                   GROUP_NAME_SPECIES,
                                   measurer.OUTPUT_NAME_EXPECTATION,
                                   recorder)
    
    def get_prereq_keys(self, var_name):
        marginal_key = (GROUP_NAME_SPECIES,
                        OUTPUT_NAME_SPARSE_MARGINAL,
                        var_name)
        return (marginal_key,)
    
    def compute(self, var_name, t_current, p_current, prereq_measurements):
        marginal_key = (GROUP_NAME_SPECIES,
                        OUTPUT_NAME_SPARSE_MARGINAL,
                        var_name)
        sparse_marginal = prereq_measurements[marginal_key]
        
        expectation = 0.0
        for (count, prob) in sparse_marginal.items():
            expectation += count*prob
        return expectation

class SpeciesStdDevMeasurer(measurer.Measurer): 
    def __init__(self, recorder):
        measurer.Measurer.__init__(self,
                                   GROUP_NAME_SPECIES,
                                   measurer.OUTPUT_NAME_STD_DEV,
                                   recorder)
    
    def get_prereq_keys(self, var_name):
        marginal_key = (GROUP_NAME_SPECIES,
                        OUTPUT_NAME_SPARSE_MARGINAL,
                        var_name)
        expectation_key = (GROUP_NAME_SPECIES,
                           measurer.OUTPUT_NAME_EXPECTATION,
                           var_name)
        return (marginal_key, expectation_key)
    
    def compute(self, var_name, t_current, p_current, prereq_measurements):
        marginal_key = (GROUP_NAME_SPECIES,
                        OUTPUT_NAME_SPARSE_MARGINAL,
                        var_name)
        sparse_marginal = prereq_measurements[marginal_key]
        
        expectation_key = (GROUP_NAME_SPECIES,
                           measurer.OUTPUT_NAME_EXPECTATION,
                           var_name)
        expectation = prereq_measurements[expectation_key]
        
        std_dev = 0.0
        for (count, prob) in sparse_marginal.items():
            std_dev += ((count - expectation) ** 2.0) * prob
        return std_dev ** 0.5
