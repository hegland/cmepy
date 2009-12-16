import unittest
from test import test_support

import numpy
from numpy.testing.utils import assert_equal

def compute_run_ends(a):
    ends = numpy.zeros(numpy.shape(a), dtype = numpy.bool)
    ends[-1] = True
    ends[:-1] = numpy.not_equal(a[:-1], a[1:])
    return ends

def sum_duplicate_keys(keys, values):
    sort_indices = numpy.argsort(keys)
    sorted_keys = keys[sort_indices]
    sorted_values = values[sort_indices]
    run_ends = compute_run_ends(sorted_keys)
    unique_keys = sorted_keys[run_ends]
    
    # since there is some finite number of unique keys we can
    # put them into the obvious bijection with the natural numbers (incl. 0) ...
    # (we need to do this to interface with numpy.bincount, since
    # if bincount is to return a length n array the values of the
    # first array must be integers from 0, .., n-1
    
    sorted_bijected_keys = numpy.zeros(numpy.shape(keys), dtype = numpy.int)
    sorted_bijected_keys[1:] = numpy.add.accumulate(run_ends[:-1])
    
    merged_values = numpy.bincount(sorted_bijected_keys,
                                   sorted_values)
    
    return (unique_keys, merged_values)
    
class KeyValueMergeTests(unittest.TestCase):
    def testSortAndSelectRepresentatives(self):
        keys = numpy.array([2, 3, 4, 2, 3, 3, 2, 4, 4])
        vals = numpy.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
        
        goal_unique_keys = numpy.array([2, 3, 4])
        goal_merged_vals = numpy.array([12, 13, 20])
        
        (unique_keys, merged_vals) = sum_duplicate_keys(keys, vals)
        
        assert_equal(unique_keys, goal_unique_keys)
        assert_equal(merged_vals, goal_merged_vals)

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(KeyValueMergeTests)
    return suite

def main():
    test_support.run_unittest(KeyValueMergeTests)

if __name__ == '__main__':
    main()