import unittest
from test import test_support

import numpy
from numpy.testing.utils import assert_almost_equal

from cmepy.statistics import Distribution
from cmepy.measurement import Measurement

class MeasurementTests(unittest.TestCase):
    def test_attribute_trickery(self):
        m = Measurement()
        m.write(0.0, Distribution({(0, ) : 1.0, (1, ) : 0.0}))
        m.write(0.5, Distribution({(0, ) : 0.5, (1, ) : 0.5}))
        m.write(1.0, Distribution({(0, ) : 0.0, (1, ) : 1.0}))
        
        assert_almost_equal(numpy.array(m.expectation),
                            [[0], [0.5], [1]])
        
        assert_almost_equal(numpy.array(m.variance),
                            [0, 0.25, 0])

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(MeasurementTests)
    return suite

def main():
    test_support.run_unittest(MeasurementTests)

if __name__ == '__main__':
    main()

        
        