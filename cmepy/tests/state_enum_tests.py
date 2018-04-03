import unittest

import numpy
from numpy.testing.utils import assert_array_equal

import cmepy.state_enum as state_enum

class StateEnumTests(unittest.TestCase):
        
    def test_simple_init(self):
        
        states = [[0, 0, 1, 1, 2, 1, 2],
                  [0, 1, 0, 1, 1, 2, 2]]
        
        enum = state_enum.create(states)
        
        enum.indices(states)
        enum.states([6, 5, 4, 3, 2, 1, 0])
        
        # verify that ordered_states are the same as manually sorting
        # with lexsort
        
        sorted_states = numpy.array(states)
        order = numpy.lexsort(sorted_states)
        sorted_states = sorted_states[:, order]
        assert_array_equal(enum.ordered_states,
                           sorted_states)
        
        # verify that the states have been enumerated with
        # respect to lexical ordering
        
        assert_array_equal(enum.states(order), states)
        
        assert_array_equal(enum.indices(states), order)
    
    def test_non_unique_indices_query(self):
        states = [[0, 0, 1, 1, 2, 7, 2],
                  [0, 1, 0, 1, 1, 0, 2]]
        
        goal_order = [0, 2, 5, 1, 3, 4, 6]
        assert_array_equal(numpy.lexsort(states),
                           goal_order)
        
        enum = state_enum.create(states)
        
        # verify that looking up the indices
        # for a collection of non-unique query states
        # returns a correctly sized index array with the correct
        # corresponding indices
        
        query_states = [[7, 2, 1, 2, 1, 7, 7, 2, 1],
                        [0, 1, 0, 2, 0, 0, 0, 2, 1]]
        
        indices = enum.indices(query_states)
        
        assert_array_equal(indices,
                           [2, 5, 1, 6, 1, 2, 2, 6, 4])
        
    def test_non_unique_contains_query(self):
        states = [[0, 0, 1, 1, 2, 7, 2],
                  [0, 1, 0, 1, 1, 0, 2]]
        
        enum = state_enum.create(states)
        
        query_states = [[-1, 7, 2, 1, 2, 9, 1, 7, 7, 2, -1, 1],
                        [-1, 0, 1, 0, 2, 9, 0, 0, 0, 2, -1, 1]]
        
        member_flags = enum.contains(query_states)
        
        goal_member_flags = [False,
                             True,
                             True,
                             True,
                             True,
                             False,
                             True,
                             True,
                             True,
                             True,
                             False,
                             True]
        
        assert_array_equal(member_flags,
                           goal_member_flags)
    
    def test_non_unique_indices_query_offset(self):
        states = [[0, 0, 1, 1, 2, 7, 2],
                  [0, 1, 0, 1, 1, 0, 2]]
        
        goal_order = [0, 2, 5, 1, 3, 4, 6]
        assert_array_equal(numpy.lexsort(states),
                           goal_order)
        
        enum = state_enum.create(states)
        
        # apply offset
        offset = 42
        enum.offset = offset
        
        # verify that looking up the indices
        # for a collection of non-unique query states
        # returns a correctly sized index array with the correct
        # corresponding indices
        
        query_states = [[7, 2, 1, 2, 1, 7, 7, 2, 1],
                        [0, 1, 0, 2, 0, 0, 0, 2, 1]]
        
        indices = enum.indices(query_states)
        
        assert_array_equal(indices - offset,
                           [2, 5, 1, 6, 1, 2, 2, 6, 4])
    
    def test_pack_unpack_distributions(self):
        states = [[0, 0, 1, 1, 2, 7, 2],
                  [0, 1, 0, 1, 1, 0, 2]]
        
        enum = state_enum.create(states)
        
        p_sparse = {(1, 1) : 0.3,
                    (7, 0) : 0.2,
                    (1, 0) : 0.5}
        
        p_dense = enum.pack_distribution(p_sparse)
        
        assert_array_equal(p_dense,
                           [0.0, 0.5, 0.2, 0.0, 0.3, 0.0, 0.0])
        
        q_sparse = enum.unpack_distribution(p_dense)
        
        for state in p_sparse:
            assert state in q_sparse
            assert p_sparse[state] == q_sparse[state]
        for state in q_sparse:
            assert state in p_sparse
            assert p_sparse[state] == q_sparse[state]
    

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(StateEnumTests)
    return suite

def main():
    unittest.run(StateEnumTests)

if __name__ == '__main__':
    main()
