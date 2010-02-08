import unittest
from test import test_support

import numpy
from numpy.testing.utils import assert_array_equal

import cmepy.new_core.lexarrayset as lexarrayset

def assert_las_equal(las1, las2):
    assert_array_equal(las1.data, las2.data)

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

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(LexArraySetTests)
    return suite

def main():
    test_support.run_unittest(LexArraySetTests)

if __name__ == '__main__':
    main()
