import unittest
from test import test_support

import itertools

import numpy
from numpy.testing.utils import assert_almost_equal

import cmepy.new_core.recorder

#
# XXX TODO ADD NORIGIN TEST CASES
# NON-NULL ZERO NORIGIN
# NULL NORIIGN
# NON-NULL NON-ZERO NORIGIN

# PRODUCT OF ALL ABOVE WITH (EXPLICIT, DEFAULT) TRANSFORMS

class RecorderTests(unittest.TestCase):
    def setUp(self):
        dims = 2
        inflated_shape = (5, )*dims
        
        def create_species_count(i):
            return lambda *x : x[i]
        
        propensities = (None, )*dims
        species_counts = [create_species_count(i) for i in xrange(dims)]
        species = ['s_%d' % i for i in xrange(dims)]
        model = {'np' : inflated_shape,
                 'propensities' : propensities,
                 'species counts' : species_counts,
                 'species' : species}
        
        goal_expected_values = {'s_0' : 100.0, 's_1' : 150.0}
        
        recorder_args = {'model' : model}
        
        target_args = {'group_name' : 'species',
                       'outputs' : ['expected value'],
                       'variable_names' : species,}
        
        dummy_p = numpy.indices(inflated_shape)[1]
        data = [(0.0, dummy_p)]
        
        def generate_cases():
            for transform in (species_counts, None):
                # XXX TODO EXTEND TO INTERESTING NON-ZERO ONES
                for norigin in (None, (0, )*dims):
                    case_recorder_args = dict(recorder_args)
                    case_recorder_args['norigin'] = norigin
                    case_target_args = dict(target_args)
                    case_target_args['transforms'] = transform
                    case_data = data
                    case_goals = {}
                    
                    for dim in xrange(dims):
                        var_name = species[dim]
                        var_goal = goal_expected_values[var_name]
                        if norigin is not None:
                            # XXX TODO FIX THIS
                            var_goal += norigin[dim]*10
                        case_goals[('species',
                                    'expected value',
                                    var_name)] = [var_goal]
                    yield (case_recorder_args,
                           case_target_args,
                           case_data,
                           case_goals)
        
        self.recorder_test_cases = generate_cases()
        
                   
    
    def testRecorder(self):
        for recorder_args, target_args, data, goals in self.recorder_test_cases:
            # test recorder init
            recorder = cmepy.new_core.recorder.CmeRecorder(**recorder_args)
            # add some targets to record
            recorder.add_target(**target_args)
            # write some dummy measurements
            for dummy_t, dummy_p in data:
                recorder.write(dummy_t, dummy_p)
            
            target_group_name = target_args['group_name']
            target_output_names = target_args['outputs']
            num_target_vars = len(target_args['variable_names'])
            target_measurement_gen = recorder.measurements(target_group_name)
            
            assert target_measurement_gen is not None
            target_measurements = [m for m in target_measurement_gen]
            assert len(target_measurements) == num_target_vars
            
            for measurement in target_measurements:
                assert measurement is not None
                assert len(measurement.times) == len(data)
                for target_output_name in target_output_names:
                    output_name = target_output_name.replace(' ', '_')
                    assert output_name in measurement.__dict__
                    output_data = measurement.__dict__[output_name]
                    assert output_data is not None
                    assert len(output_data) == len(data)
                    for i, measured_t in enumerate(measurement.times):
                        assert measured_t == data[i][0]
                        goal_data = goals[(target_group_name,
                                          target_output_name,
                                          measurement.name)][i]
                        assert_almost_equal(output_data[i],
                                            goal_data)
    
        
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(RecorderTests)
    return suite

def main():
    test_support.run_unittest(RecorderTests)

if __name__ == '__main__':
    main()
