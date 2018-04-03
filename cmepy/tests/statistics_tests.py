import unittest
import numpy
from numpy.testing.utils import assert_almost_equal


from cmepy import statistics

class StatisticsTests(unittest.TestCase):
    def test_one_dee_distributions(self):
        p_1 = {(1, ) : 0.1,
               (2, ) : 0.2,
               (3, ) : 0.3,
               (4, ) : 0.4,}
        
        f = lambda state : (0, )
        p = statistics.map_distribution(f, p_1)
        assert len(p) == 1
        assert ((0, ) in p) and (p[(0, )] == 1.0)
        
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
        
        f = lambda state : (state[0], )
        p = statistics.map_distribution(f, p_2)
        assert len(p) == 3
        assert ((0, ) in p) and (p[(0, )] == 0.5)
        assert ((1, ) in p) and (p[(1, )] == 0.2)
        assert ((3, ) in p) and (p[(3, )] == 0.3)
        
        mu = statistics.expectation(p_2)
        mu_goal = numpy.asarray([1.1, 1.2])
        assert_almost_equal(mu, mu_goal)
        
        cov = statistics.covariance(p_2)
        cov_goal = -1.1*-1.2*0.2 -1.1*-0.2*0.3 -0.1*-1.2*0.2 +1.9*1.8*0.3
        assert_almost_equal(cov, cov_goal)
    
    def test_dense_conversions(self):
        p = statistics.Distribution()
        
        p_dense = numpy.array(
            [
                [0.0, 0.2, 0.1, 0.0],
                [0.3, 0.0, 0.1, 0.3],
            ]
        )
        
        p.from_dense(p_dense)
        
        assert len(p) == 5
        
        assert p[(0, 1)] == 0.2
        assert p[(0, 2)] == 0.1
        
        assert p[(1, 0)] == 0.3
        assert p[(1, 2)] == 0.1
        assert p[(1, 3)] == 0.3
        
        p.from_dense(p_dense, origin = (-1, 3))
        
        assert len(p) == 5
        
        assert p[(-1, 4)] == 0.2
        assert p[(-1, 5)] == 0.1
        
        assert p[(0, 3)] == 0.3
        assert p[(0, 5)] == 0.1
        assert p[(0, 6)] == 0.3
    
    def test_distribution_addition(self):
        a = statistics.Distribution(
            {
                (0, 0) : 0.4,
                (0, 1) : 0.3,
                (1, 0) : 0.2,
                (1, 1) : 0.1
            }
        )
        
        b = statistics.Distribution(
            {
                (0, -1) : 0.1,
                (0, 0) : 0.2,
                (0, 1) : 0.3,
                (0, 2) : 0.4
            }
        )
        
        c = a + b
        
        assert len(c) == 6
        assert_almost_equal(c[(0, -1)], 0.1)
        assert_almost_equal(c[(0, 0)], 0.6)
        assert_almost_equal(c[(0, 1)], 0.6)
        assert_almost_equal(c[(1, 0)], 0.2)
        assert_almost_equal(c[(1, 1)], 0.1)
        assert_almost_equal(c[(0, 2)], 0.4)
    
    def test_distribution_subtraction(self):
        a = statistics.Distribution(
            {
                (0, 0) : 0.4,
                (0, 1) : 0.3,
                (1, 0) : 0.2,
                (1, 1) : 0.1
            }
        )
        
        b = statistics.Distribution(
            {
                (0, -1) : 0.1,
                (0, 0) : 0.2,
                (0, 1) : 0.3,
                (0, 2) : 0.4
            }
        )
        
        c = b - a
        
        assert len(c) == 6
        assert_almost_equal(c[(0, -1)], 0.1)
        assert_almost_equal(c[(0, 0)], -0.2)
        assert_almost_equal(c[(0, 1)], 0.0)
        assert_almost_equal(c[(1, 0)], -0.2)
        assert_almost_equal(c[(1, 1)], -0.1)
        assert_almost_equal(c[(0, 2)], 0.4)
    
    def test_distribution_scalar_multiplication(self):
        a = statistics.Distribution(
            {
                (0, 0) : 0.4,
                (0, 1) : 0.3,
                (1, 0) : 0.2,
                (1, 1) : 0.1
            }
        )
        
        b = - 3.2
        
        # left mult with scalar
        c = a * b
        
        assert len(c) == 4
        assert_almost_equal(c[(0, 0)], 0.4*b)
        assert_almost_equal(c[(0, 1)], 0.3*b)
        assert_almost_equal(c[(1, 0)], 0.2*b)
        assert_almost_equal(c[(1, 1)], 0.1*b)
        
        # right mult with scalar
        c = b * a
        
        assert len(c) == 4
        assert_almost_equal(c[(0, 0)], 0.4*b)
        assert_almost_equal(c[(0, 1)], 0.3*b)
        assert_almost_equal(c[(1, 0)], 0.2*b)
        assert_almost_equal(c[(1, 1)], 0.1*b)
    
    def test_distribution_convex_combination(self):
        a = statistics.Distribution(
            {
                (0, 0) : 0.4,
                (0, 1) : 0.3,
                (1, 0) : 0.2,
                (1, 1) : 0.1
            }
        )
        
        b = statistics.Distribution(
            {
                (0, -1) : 0.1,
                (0, 0) : 0.2,
                (0, 1) : 0.3,
                (0, 2) : 0.4
            }
        )
        
        tau = 0.73111
        
        # let c be a convex combination of a and b, written in a wierd way.
        # this should test scalar mult, addition, unary negation operator
        # etc.
        c = a * tau + b + tau * (-b)
        
        assert len(c) == 6
        assert_almost_equal(c[(0, -1)], 0.1 * (1.0 - tau))
        assert_almost_equal(c[(0, 0)], 0.4 * tau + 0.2 * (1.0 - tau))
        assert_almost_equal(c[(0, 1)], 0.3 * tau + 0.3 * (1.0 - tau))
        assert_almost_equal(c[(1, 0)], 0.2 * tau)
        assert_almost_equal(c[(1, 1)], 0.1 * tau)
        assert_almost_equal(c[(0, 2)], 0.4 * (1.0 - tau))
    
    def test_norms_and_metrics(self):
        a = statistics.Distribution(
            {
                (0, 0) : 0.4,
                (0, 1) : 0.3,
                (1, 0) : 0.2,
                (1, 1) : 0.1
            }
        )
        
        for p in numpy.linspace(0.1, 2.5, 25):
            assert_almost_equal(a.lp_distance(a, p), 0.0)
            assert a.lp_norm(p) > 0.0
        
        assert numpy.isinf(a.kl_divergence(statistics.Distribution()))
        assert_almost_equal(a.kl_divergence(a), 0.0)
        
        
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(StatisticsTests)
    return suite

def main():
    unittest.run(StatisticsTests)

if __name__ == '__main__':
    main()
