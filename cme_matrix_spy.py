import numpy
import scipy.sparse
import pylab

import perturb
import rank_approx

import cmepy.new_core.cme_solver as cme_solver
import cmepy.core.matrix_cme as matrix_cme

def transform_coo_matrix_basis(coo_matrix, new_indices):
    n, n = coo_matrix.shape
    assert numpy.size(new_indices) == n
    new_matrix_data = (coo_matrix.data,
                       (new_indices[coo_matrix.row],
                        new_indices[coo_matrix.col]))
    return scipy.sparse.coo_matrix(new_matrix_data,
                                   coo_matrix.shape)

def create_nice_matrix(model, f=None):
    flux_data = cme_solver.create_flux_data(model)
    matrix = matrix_cme.gen_sparse_matrix(model['np'], flux_data)
    if f is not None:
        # change basis so that states are enumerated in order of increasing
        # value under the transform f
        indices = numpy.arange((numpy.product(model['np'],)))
        states = [numpy.ravel(i) for i in numpy.indices(model['np'])]
        f_states = f(*states)
        new_order = numpy.argsort(f_states)
        new_indices = indices[new_order]
                
        matrix = transform_coo_matrix_basis(matrix, new_indices)
        
    return matrix

def spy_sparse_matrix(sparse_matrix):
    matrix_dense = sparse_matrix.todense()
    
    pylab.figure()
    pylab.spy(matrix_dense)

def nice_spy(model, f):
    spy_sparse_matrix(create_nice_matrix(model, f))

def main():
    epsilon = 0.1
    s_0 = 10
    e_0 = 2
    model, slow_r, fast_r = rank_approx.get_catalyser_model(s_0,
                                                            e_0,
                                                            epsilon)
    (slow_model, fast_model) = perturb.decompose_model(model,
                                                       (slow_r, fast_r))
    #
    # second function here is a horrible, dirty hack.
    # what we actually mean is 'sort by x[0] then by x[2] then by x[1]
    # equivalently.... we could reorder the x's
    # try 0, 2, 1 or 1, 2, 0 .....
    # this is to gain a block diagonal representation
    
    f_range = [None, lambda *x : x[0]*1000 + x[2]*10 + x[1]]
    for f in f_range:
        nice_spy(slow_model, f)
        pylab.title('slow! f = %s' % str(f))
        nice_spy(fast_model, f)
        pylab.title('fast! f = %s' % str(f))
    
    pylab.show()

if __name__ == '__main__':
    main()
