"""
hook allowing all unit tests to be run via setuptools
"""

import unittest


# root package to gather unit tests from
ROOT_PACKAGE = 'cmepy'

# package structure, leaves must contain suite() method that returns
# the test suite for the sub package
SUB_PACKAGES = {
    ROOT_PACKAGE : [
        'tests',
    ],
    'tests' : [
        'ode_solver_tests',
        'solver_tests',
        'recorder_tests',
        'domain_tests',
        'state_enum_tests',
        'lexarrayset_tests',
        'statistics_tests',
        'measurement_tests',
        'model_tests',
    ],
}



def import_sub_package(sub_package_chain):
    """
    import_sub_package(['A', 'B', 'C']) -> module
    
    where module is the sub package A.B.C
    """
    
    name = '.'.join(sub_package_chain)
    fromlist = [sub_package_chain[-1]]
    return __import__(name, fromlist = fromlist)

def gather_suites(sub_packages, root_package):
    """
    gather_suites(sub_packages, root_package) -> set of unittest.TestSuite
    """
    test_suites = set()
    
    def dfs_add_tests(chain):
        """
        dfs_add_tests(chain)
        
        recursively explore package structure using dfs, loading
        leaf sub packages, and adding the results of their suite()
        method to the set test_suites.
        """
        head = chain[-1]
        if head not in sub_packages:
            sub_package = import_sub_package(chain)
            print '\t+ %s' % str(sub_package)
            try:
                suite = sub_package.suite()
                test_suites.add(suite)
            except AttributeError:
                detail = 'sub package \'%s\' has no \'suite()\' method' % head
                print '\t  -- WARNING : %s, ignoring' % detail
        else:
            for sub_package in sub_packages[head]:
                dfs_add_tests(chain + [sub_package])
    
    dfs_add_tests([root_package])
    
    return test_suites

def additional_tests(sub_packages, root_package):
    """
    additional_tests() -> unittest.TestSuite
    """
    
    print ''
    print '-- gathering test suites :'
    print ''
    test_suite = gather_suites(sub_packages, root_package)
    print ''
    print '-- running test suites :'
    print ''
    all_test_suite = unittest.TestSuite(test_suite)
    
    return all_test_suite

def main():
    """
    gathers and runs all the tests
    """
    suite = additional_tests(SUB_PACKAGES, ROOT_PACKAGE)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
if __name__ == '__main__':
    main()
