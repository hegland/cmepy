import unittest

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
        accumulator = block_diagonal.Accumulator()
        k = 0
        while k<n:
            block_size = numpy.random.randint(m)+1
            block_size = min(block_size, n-k)
            block = numpy.random.random((block_size, block_size))+1
            accumulator.add_dense_block(block,
                                        k,
                                        k,
                                        block_size)
            k += block_size
        yield accumulator.to_coo_matrix()

class BlockDiagonalTests(unittest.TestCase):
    def setUp(self):
        self.n = 20
        m = 5
        runs = 30
        #self.matrices = [scipy.sparse.coo_matrix(numpy.array([[0.0]]))]
        self.matrices = []
        for run, matrix in enumerate(random_block_diag_matrices(self.n, m)):
            if run >= runs:
                break
            self.matrices.append(matrix)
    
    def testBlockDiagExpm(self):
        for matrix in self.matrices:
            dense_matrix = matrix.todense()
            t = 20.0
            goal_expm = scipy.linalg.expm(dense_matrix*t)
            
            block_diag = block_diagonal.from_sparse_matrix(matrix)
            block_diag_expm = block_diagonal.expm(block_diag, t)
            sparse_expm = block_diagonal.to_sparse(block_diag_expm)
            test_expm = sparse_expm.todense()
            
            #error = scipy.linalg.norm(test_expm - goal_expm)
            #print 'norm of error : '+str(error)
            #rel_error = error / scipy.linalg.norm(goal_expm)
            #print 'norm of error (relative to goal norm) : '+str(rel_error)
            
            error = numpy.ravel(test_expm - goal_expm)
            tol = 1.e-15
            bad_elements = numpy.compress(numpy.abs(error)>=tol, error)
            print 'bad_elements : '
            print str(bad_elements)
            try:
                numpy.testing.utils.assert_almost_equal(test_expm, goal_expm)
            except AssertionError:
                import pylab
                tol = 1.0e-10
                pylab.spy(numpy.abs(test_expm - goal_expm)>tol)
                pylab.show()
                break
    
    def testBlockDiagSVD(self):
        for matrix in self.matrices:
            # traditional dense routine
            dense_matrix = matrix.todense()
            k = numpy.random.randint(1, self.n+1)
            print 'using rank k = '+str(k)
            u, s, vh = scipy.linalg.svd(dense_matrix)
            f_goal = u[:, :k]
            e_goal = numpy.dot(numpy.diag(s[:k]), vh[:k, :])
            
            # new block diagonal routine
            block_diag = block_diagonal.from_sparse_matrix(matrix)
            block_diag_svd = block_diagonal.block_svd(block_diag)
            test_approx = block_diagonal.to_sparse_rank_k_approx(block_diag_svd,
                                                                 k)
            e_test, f_test = test_approx
            numpy.testing.utils.assert_almost_equal(e_test.todense(), e_goal)
            numpy.testing.utils.assert_almost_equal(f_test.todense(), f_goal)
if __name__ == '__main__':
    unittest.main()