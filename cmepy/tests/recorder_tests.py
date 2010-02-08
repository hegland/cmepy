import unittest
from test import test_support

import itertools

import numpy
from numpy.testing.utils import assert_almost_equal

import cmepy.recorder

class RecorderTests(unittest.TestCase):
    # XXX TODO
    pass    
        
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(RecorderTests)
    return suite

def main():
    test_support.run_unittest(RecorderTests)

if __name__ == '__main__':
    main()
