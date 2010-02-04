import unittest
from test import test_support

import itertools

import numpy
from numpy.testing.utils import assert_almost_equal

import cmepy.recorder

class RecorderTests(unittest.TestCase):
    def test_recorder_a(self):
        """
        test some statistical operations
        """
        p_1 = {(1, ) : 0.1,
               (2, ) : 0.2,
               (3, ) : 0.3,
               (4, ) : 0.4,}
        
        rec = cmepy.recorder.create((('foo',),))
        
        rec.write(1.0, p_1)
        
        ev = rec['foo'].expected_value
        assert_almost_equal(ev, numpy.array([[3.0]]))
        
        var = rec['foo'].variance
        assert_almost_equal(var, numpy.array([1.0]))
        std_dev = rec['foo'].standard_deviation
        assert_almost_equal(std_dev, numpy.array([1.0]))
        
        
    def test_recorder_b(self):
        """
        test out some more statistical operations
        """
        p_1 = {(1, ) : 0.1,
               (2, ) : 0.2,
               (3, ) : 0.3,
               (4, ) : 0.4,}
        
        # map states via pairwise aggregation
        even = lambda *x : x[0]/2
        odd = lambda *x : x[0]/2 + 1
        
        rec = cmepy.recorder.create((('even', 'odd'), (even, odd)))
        
        rec.write(1.0, p_1)
        
        even_ev = rec['even'].expected_value
        assert_almost_equal(even_ev, numpy.array([[1.3]]))
        odd_ev = rec['odd'].expected_value
        assert_almost_equal(odd_ev, numpy.array([[2.3]]))
        
        cov = rec[('even', 'odd')].covariance
        assert_almost_equal(cov, numpy.array([0.41]))
        
        
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(RecorderTests)
    return suite

def main():
    test_support.run_unittest(RecorderTests)

if __name__ == '__main__':
    main()
