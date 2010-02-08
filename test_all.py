"""
hook allowing all unit tests to be run via setuptools
"""

import unittest

def additional_tests():
    test_suites = set()
    
    from cmepy.core import cme_tests
    test_suites.add(cme_tests.suite())
    
    from cmepy.new_core import ode_solver_tests
    test_suites.add(ode_solver_tests.suite())
    
    from cmepy.new_core import cme_solver_tests
    test_suites.add(cme_solver_tests.suite())
    
    from cmepy.new_core import recorder_tests
    test_suites.add(recorder_tests.suite())
    
    from cmepy.new_core import domain_tests
    test_suites.add(domain_tests.suite())
    
    from cmepy.new_core import state_enum_tests
    test_suites.add(state_enum_tests.suite())
    
    all_test_suite = unittest.TestSuite(test_suites)
    
    return all_test_suite

if __name__ == '__main__':
    suite = additional_tests()
    unittest.TextTestRunner(verbosity=2).run(suite)


