import unittest

import itertools

import numpy
import numpy.random
import numpy.testing.utils
import scipy.sparse
import scipy.linalg


import block_diagonal

def random_block_diag_matrices(n, m):
    """
    returns a generator for an infinite sequence of
    nxn sparse coo matrices with 1 <= block_size <= m for each block
    """
    while True:
        accumulator = block_diagonal.Accumulator((n, n))
        k = 0
        while k<n:
            block_size = numpy.random.randint(m)+1
            block_size = min(block_size, n-k)
            block = numpy.random.random((block_size, block_size))
            col_means = numpy.add.reduce(block, axis=1) / block_size
            col_means = numpy.reshape(col_means, (1, block_size))
            block = block - col_means
            accumulator.add_dense_block(block,
                                        k,
                                        k)
            k += block_size
        yield accumulator.to_coo_matrix()

class BlockDiagonalTests(unittest.TestCase):
    def setUp(self):
        self.n = 400
        m = 20
        runs = 30
        #self.matrices = [scipy.sparse.coo_matrix(numpy.array([[0.0]]))]
        self.matrices = []
        for run, matrix in enumerate(random_block_diag_matrices(self.n, m)):
            self.matrices.append(matrix)
            if run+1 == runs:
                break
    
    def testBlockDiagExpm(self):
        for matrix in self.matrices:
            dense_matrix = matrix.todense()
            t = 1.0
            goal_expm = scipy.linalg.expm(dense_matrix*t)
            
            block_diag = block_diagonal.from_sparse_matrix(matrix)
            block_diag_expm = block_diagonal.expm(block_diag, t)
            sparse_expm = block_diagonal.to_sparse(block_diag_expm)
            test_expm = sparse_expm.todense()
            
            numpy.testing.utils.assert_almost_equal(test_expm, goal_expm)
            
            """
            try:
                numpy.testing.utils.assert_almost_equal(test_expm, goal_expm)
            except AssertionError:
                print 'numpy assert almost equals failure'
                print 'dumping detailed information : '
                tol = 1.0e-7
                print 'test shape : %s' % str(numpy.shape(test_expm))
                print 'goal shape : %s' % str(numpy.shape(goal_expm))
                for x, y in itertools.izip(numpy.ravel(test_expm),
                                           numpy.ravel(goal_expm)):
                    print ('\t%g ~ %g\t' % (x, y)),
                    if abs(x-y) > tol:
                        print '***FAIL***'
                    else:
                        print 'OK'
                #import pylab
                #tol = 1.0e-10
                #pylab.spy(numpy.abs(test_expm - goal_expm)>tol)
                #pylab.show()
                #break
                print 'continuing ...'
            """
    
    def testBlockDiagSVD(self):
        compare_with_dense = False
        for matrix in self.matrices:
            k = numpy.random.randint(1, self.n+1)
            
            # new block diagonal routine
            block_diag = block_diagonal.from_sparse_matrix(matrix)
            block_diag_svd = block_diagonal.block_svd(block_diag)
            test_approx = block_diagonal.to_sparse_rank_k_approx(block_diag_svd,
                                                                 k)
            e_test, f_test = test_approx
            if compare_with_dense:
                print 'comparing with dense for k = %d' % k
                # traditional dense routine
                dense_matrix = matrix.todense()
                u, s, vh = scipy.linalg.svd(dense_matrix)
                f_goal = u[:, :k]
                e_goal = numpy.dot(numpy.diag(s[:k]), vh[:k, :])
                test_approx = numpy.dot(f_test.todense(), e_test.todense())
                goal_approx = numpy.dot(f_goal, e_goal)
                numpy.testing.utils.assert_almost_equal(test_approx, goal_approx)

def main():
    unittest.main()

if __name__ == '__main__':
    import cProfile, pstats
    PROFILE_FILE = 'block_diagonal_tests.profile'
    cProfile.run('main()', PROFILE_FILE)
    stats = pstats.Stats(PROFILE_FILE)
    stats.sort_stats('cumulative').print_stats(30)