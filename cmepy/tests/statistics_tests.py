import unittest
from test import test_support
import numpy
from numpy.testing.utils import assert_almost_equal


from cmepy import statistics

class StatisticsTests(unittest.TestCase):
    def test_one_dee_distributions(self):
        p_1 = {(1, ) : 0.1,
               (2, ) : 0.2,
               (3, ) : 0.3,
               (4, ) : 0.4,}
        
        f = lambda state : 0
        p = statistics.map_distribution(f, p_1)
        assert len(p) == 1
        assert (0 in p) and (p[0] == 1.0)
        
        mu = statistics.expectation(p_1)
        mu_goal = numpy.array([3.0])
        assert_almost_equal(mu, mu_goal)
        
        sigma_squared = statistics.variance(p_1)
        sigma_squared_goal = numpy.asarray([1.0])
        assert_almost_equal(sigma_squared, sigma_squared_goal)
    
    def test_two_dee_distributions(self):
        p_2 = {(0, 0) : 0.2,
               (0, 1) : 0.3,
               (1, 0) : 0.2,
               (3, 3) : 0.3,}
        
        f = lambda state : state[0]
        p = statistics.map_distribution(f, p_2)
        assert len(p) == 3
        assert (0 in p) and (p[0] == 0.5)
        assert (1 in p) and (p[1] == 0.2)
        assert (3 in p) and (p[3] == 0.3)
        
        mu = statistics.expectation(p_2)
        mu_goal = numpy.asarray([1.1, 1.2])
        assert_almost_equal(mu, mu_goal)
        
        cov = statistics.covariance(p_2)
        cov_goal = -1.1*-1.2*0.2 -1.1*-0.2*0.3 -0.1*-1.2*0.2 +1.9*1.8*0.3
        assert_almost_equal(cov, cov_goal)
        
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(StatisticsTests)
    return suite

def main():
    test_support.run_unittest(StatisticsTests)

if __name__ == '__main__':
    main()
