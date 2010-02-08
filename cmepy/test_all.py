"""
hook allowing all unit tests to be run via setuptools
"""

import unittest

def additional_tests():
    
    from cmepy.core import cme_tests
    cme_test_suite = cme_tests.suite()
    
    from cmepy.solver.core_cme_solver import cme_solver_tests
    cme_solver_test_suite = cme_solver_tests.suite()
    
    all_test_suite = unittest.TestSuite([cme_test_suite,
                                         cme_solver_test_suite,])
    
    return all_test_suite

if __name__ == '__main__':
    suite = additional_tests()
    unittest.TextTestRunner(verbosity=2).run(suite)


