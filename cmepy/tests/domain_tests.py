import unittest

import itertools

import numpy
from numpy.testing.utils import assert_array_equal

import cmepy.domain as domain


class DomainTests(unittest.TestCase):
    def test_rect_domain(self):
        shape = (3, 1, 3, 8)
        states = domain.from_rect(shape)
        goal_states = numpy.reshape(numpy.indices(shape), (len(shape), -1))
        assert_array_equal(states, goal_states)
    
    def test_rect_offset_domain(self):
        shape = (4, 7, 3)
        origin = (3, -2, 2)
        states = domain.from_rect(shape, origin=origin)
        vect_origin = numpy.asarray(origin)[:, numpy.newaxis]
        goal_states = numpy.reshape(numpy.indices(shape), (len(shape), -1))
        goal_states += vect_origin
        assert_array_equal(states, goal_states)
    
    def test_sparse_domain_from_dict(self):
        p_0 = {(4, 3) : 0.1,
               (9, 4) : 0.3,
               (11, -1) : 0.2,
               (44, 44) : 0.4}
        
        states = domain.from_iter(p_0)
        
        goal_states = list(p_0)
        for i, goal_state in enumerate(goal_states):
            assert_array_equal(states[:, i], goal_state)
    
    def test_sparse_domain_from_set(self):
        sparse_states = set([(1, 2, 3),
                             (2, 3, 4),
                             (3, 4, 5),
                             (4, 5, 6),
                             (2, 3, 4),
                             (2, 3, 4),
                             (1, 2, 3)])
        
        states = domain.from_iter(sparse_states)
        
        goal_states = list(sparse_states)
        for i, goal_state in enumerate(goal_states):
            assert_array_equal(states[:, i], goal_state)
    
    def test_iter_from_dense_states(self):
        sparse_states = set([(1, 2, 3),
                             (2, 3, 4),
                             (3, 4, 5),
                             (4, 5, 6),
                             (2, 3, 4),
                             (2, 3, 4),
                             (1, 2, 3)])
        
        dense_states = domain.from_iter(sparse_states)
        state_iter = domain.to_iter(dense_states)
        for (state, goal_state) in itertools.izip(state_iter, sparse_states):
            assert state == goal_state
    
    def test_sparse_domain_from_mapping(self):
        p_0 = {(4, 3) : 0.1,
               (9, 4) : 0.3,
               (11, -1) : 0.2,
               (44, 44) : 0.4}
        
        states, values = domain.from_mapping(p_0)
        
        goal_states = list(p_0)
        for i, goal_state in enumerate(goal_states):
            assert_array_equal(states[:, i], goal_state)
            assert p_0[goal_state] == values[i]
    
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(DomainTests)
    return suite

def main():
    unittest.run(DomainTests)

if __name__ == '__main__':
    main()
