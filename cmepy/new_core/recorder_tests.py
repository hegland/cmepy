import unittest
from test import test_support

import itertools

import numpy
from numpy.testing.utils import assert_almost_equal

import cmepy.new_core.recorder

class RecorderTests(unittest.TestCase):
    def testSpeciesCountExpectation(self):
        dims = 2
        inflated_shape = (5, )*dims
        
        def create_species_count(i):
            return lambda *x : x[i]
        
        propensities = (None, )*dims
        species_counts = [create_species_count(i) for i in xrange(dims)]
        species = ['s_%d' % i for i in xrange(dims)]
        model = {'propensities' : propensities,
                 'species counts' : species_counts,
                 'species' : species}
        
        recorder = cmepy.new_core.recorder.CmeRecorder(model)
        
        recorder.add_target('species',
                            outputs = ['expected value', 'marginal'],
                            variable_names = model['species'],
                            transforms = model['species counts'])
        
        dummy_p = numpy.indices(inflated_shape)[1]
        
        recorder.write(0.0, dummy_p)
        
        species_measurements = [m for m in recorder.measurements('species')]
        assert species_measurements is not None
        assert len(species_measurements) == dims
        
        goals = {'s_0' : 100.0, 's_1' : 150.0}
        for measurement in species_measurements:
            assert measurement is not None
            assert len(measurement.times) == 1
            assert measurement.times[0] == 0.0
            assert len(measurement.expected_value) == 1            
            goal_expectation = goals[measurement.name]
            assert_almost_equal(measurement.expected_value[0],
                                goal_expectation)
        
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(RecorderTests)
    return suite

def main():
    test_support.run_unittest(RecorderTests)

if __name__ == '__main__':
    main()
