import unittest
import numpy.testing.utils

import lexarrayset
import fsp_structures

def assert_las_equal(las1, las2):
    assert_array_equal(las1.data, las2.data)

def assert_array_equal(states1, states2):
    numpy.testing.utils.assert_array_equal(states1, states2)

def verify_state_enum(state_enum,
                      goal_unordered_states,
                      goal_ordered_states,
                      goal_indices):
    assert_array_equal(state_enum.unordered_states,
                       goal_unordered_states)
    assert_array_equal(state_enum.ordered_states,
                       goal_ordered_states)
    assert_array_equal(state_enum.index,
                       goal_indices)

class LexArraySetTests(unittest.TestCase):
    def testEmptyLas(self):
        """
        ensure the empty LexArraySet behaves in some reasonable fashion
        """
        empty = lexarrayset.empty(3)
        las1 = lexarrayset.create([[0, 0], [1, 1], [2, 3]])
        union1 = empty.union(las1)
        union2 = las1.union(empty)
        union3 = empty.union(empty)
        assert_las_equal(union1, las1)
        assert_las_equal(union2, las1)
        assert_las_equal(union3, empty)
        
        intersection1 = empty.intersection(las1)
        intersection2 = las1.intersection(empty)
        intersection3 = empty.intersection(empty)
        assert_las_equal(intersection1, empty)
        assert_las_equal(intersection2, empty)
        assert_las_equal(intersection3, empty)
        
        difference1 = empty.difference(las1)
        difference2 = las1.difference(empty)
        assert_las_equal(difference1, empty)
        assert_las_equal(difference2, las1)
        
        split1_intersect, split1_diff = empty.split(las1)
        assert_las_equal(split1_intersect, empty)
        assert_las_equal(split1_diff, empty)
        split2_intersect, split2_diff = las1.split(empty)
        assert_las_equal(split2_intersect, empty)
        assert_las_equal(split2_diff, las1)
        split3_intersect, split3_diff = empty.split(empty)
        assert_las_equal(split3_intersect, empty)
        assert_las_equal(split3_diff, empty)
        
    def testNonUniqueArgMember(self):
        """
        test how the non-unique membership operation performs
        """
        arr1 = numpy.array([[1, 1, 1], [1, 1, -1], [0, 0, 0]])
        las2 = lexarrayset.create([[0, 1], [0, 1], [0, 0]])
        members = lexarrayset.nonunique_member(arr1, las2.data)
        assert_array_equal(members, [True,
                                     True,
                                     False])
        
        arr1 = numpy.array([[1, 3, 3, 1, 1, 1, 1],
                            [2, 2, 2, 2, 1, 1, -1],
                            [3, 1, 1, 3, 0, 0, 0]])
        las2 = lexarrayset.create([[0, 1, -1, 1],
                                   [0, 1, -1, 2],
                                   [0, 0, -1, 3]])
        members = lexarrayset.nonunique_member(arr1, las2.data)
        assert_array_equal(members, [True,
                                     False,
                                     False,
                                     True,
                                     True,
                                     True,
                                     False, ])

class StateEnumerationTests(unittest.TestCase):
    def testEmptyInitialisation(self):
        # create an empty set of 7-tuples of integers
        initial_states = lexarrayset.empty(7)
        state_enum = fsp_structures.StateEnumeration(dim=7)
        
        # extend it by some non-empty set
        sigma = lexarrayset.create([[1], [2], [3], [4], [5], [6], [7]])
        state_enum.extend(sigma)
        verify_state_enum(state_enum,
                          numpy.array([[1], [2], [3], [4], [5], [6], [7]]),
                          numpy.array([[1], [2], [3], [4], [5], [6], [7]]),
                          numpy.array([0]))
        
        # then extend it by an empty set
        state_enum.extend(lexarrayset.empty(7))
        verify_state_enum(state_enum,
                          numpy.array([[1], [2], [3], [4], [5], [6], [7]]),
                          numpy.array([[1], [2], [3], [4], [5], [6], [7]]),
                          numpy.array([0]))
        
    def testExtension(self):
        
        initial_states = lexarrayset.create([[0], [0], [0]])
        state_enum = fsp_structures.StateEnumeration(dim = 3)
        state_enum.extend(initial_states)
        
        verify_state_enum(state_enum,
                          numpy.array([[0], [0], [0]]),
                          numpy.array([[0], [0], [0]]),
                          numpy.array([0]))
        
        sigma = lexarrayset.create([[0], [1], [0]])
        state_enum.extend(sigma)
        
        verify_state_enum(state_enum,
                          numpy.array([[0, 0], [0, 1], [0, 0]]),
                          numpy.array([[0, 0], [0, 1], [0, 0]]),
                          numpy.array([0, 1]))
        
        # nb las implicitly sorts its arguments, so this makes figuring out the
        # goal unordered_state member of state_enum a touch tricky if the
        #argument used to create sigma isn't already sorted
        sigma = lexarrayset.create([[1, 0], [0, 1], [0, 1]])
        state_enum.extend(sigma)
        verify_state_enum(state_enum,
                          numpy.array([[0, 0, 1, 0], [0, 1, 0, 1], [0, 0, 0, 1]]),
                          numpy.array([[0, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 1]]),
                          numpy.array([0, 2, 1, 3]))

class ShiftStateTrackerTests(unittest.TestCase):
    def testExtension(self):
        
        interior = lexarrayset.create([[0, 1], [0, 0]])
        boundary = lexarrayset.create([[2], [0]])
        shift = numpy.asarray([1, 0])
        
        sst = fsp_structures.ShiftStateTracker(shift, dim=2)
        # artificially modify states ...
        sst.interior = interior
        sst.boundary = boundary
        
        goal_interior = lexarrayset.create([[0, 1], [0, 0]])
        goal_boundary = lexarrayset.create([[2], [0]])
        assert_las_equal(sst.interior, goal_interior)
        assert_las_equal(sst.boundary, goal_boundary)
        assert sst.interior.intersection(sst.boundary).size == 0
        
        sigma = lexarrayset.create([[3], [0]])
        sst.extend(sigma)
        
        goal_interior = lexarrayset.create([[0, 1, 2], [0, 0, 0]])
        goal_boundary = lexarrayset.create([[3], [0]])
        assert_las_equal(sst.interior, goal_interior)
        assert_las_equal(sst.boundary, goal_boundary)
        assert sst.interior.intersection(sst.boundary).size == 0
        
        sigma = lexarrayset.create([[0], [1]])
        sst.extend(sigma)
        
        goal_interior = lexarrayset.create([[0, 1, 2], [0, 0, 0]])
        goal_boundary = lexarrayset.create([[3, 0], [0, 1]])
        assert_las_equal(sst.interior, goal_interior)
        assert_las_equal(sst.boundary, goal_boundary)
        assert sst.interior.intersection(sst.boundary).size == 0
        
        sigma = lexarrayset.create([[1], [1]])
        sst.extend(sigma)
        
        goal_interior = lexarrayset.create([[0, 1, 2, 0], [0, 0, 0, 1]])
        goal_boundary = lexarrayset.create([[3, 1], [0, 1]])
        assert_las_equal(sst.interior, goal_interior)
        assert_las_equal(sst.boundary, goal_boundary)
        assert sst.interior.intersection(sst.boundary).size == 0
        

if __name__ == '__main__':
    unittest.main()
