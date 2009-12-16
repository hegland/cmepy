import unittest
from test import test_support

import numpy
import numpy.random
from numpy.testing.utils import assert_equal

def compute_run_ends(a):
    """
    maps flat array a to a boolean array e, where
    e[i] is defined by (a[i] != a[i+1])
    
    the last element of e is defined to be True
    """
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

def create_random_data(num_entries, num_keys):
    min_key = -100000
    max_key = 100000
    key_domain = numpy.random.randint(min_key, max_key, size = (num_keys, ))
    
    keys = key_domain[numpy.random.randint(0, num_keys, size = (num_entries, ))]
    values = numpy.random.random(size = (num_entries,))
    return (keys, values)

def dict_sum_duplicate_keys(keys, values):
    import itertools
    d = {}
    for k, v in itertools.izip(keys, values):
        d[k] = d.get(k, 0.0)+v
    return d

def compare():
    num_runs = 100
    for run in xrange(num_runs):
        print 'run %d of %d' % (run+1, num_runs)
        keys, values = create_random_data(100000, 100)
        result_dict = dict_sum_duplicate_keys(keys, values)
        result_vect = sum_duplicate_keys(keys, values)
    
def main():
    test_support.run_unittest(KeyValueMergeTests)

if __name__ == '__main__':
    #main()
    import cProfile, pstats
    PROFILE_FILE = 'vector_key_value_merge.profile'
    cProfile.run('compare()', PROFILE_FILE)
    stats = pstats.Stats(PROFILE_FILE)
    stats.sort_stats('cumulative').print_stats(30)
    