def create_change_of_basis_matrices(permutation):
    n = numpy.size(permutation)
    shape = (n, )*2
    forward_row = numpy.array(permutation)
    forward_col = numpy.arange(n)
    forward_val = numpy.ones((n, ))
    forward = scipy.sparse.coo_matrix((forward_val,
                                      (forward_row,
                                       forward_col)),
                                       shape).tocsr()
    inverse_row = numpy.arange(n)
    inverse_col = numpy.array(permutation)
    inverse_val = numpy.ones((n, ))
    inverse = scipy.sparse.coo_matrix((inverse_val,
                                      (inverse_row,
                                       inverse_col)),
                                       shape).tocsr()
    return (forward, inverse)

def create_ordered_basis(model, f):
    # change basis so that states are enumerated in order of increasing
    # value under the transform f
    indices = numpy.arange((numpy.product(model['np'],)))
    states = [numpy.ravel(i) for i in numpy.indices(model['np'])]
    f_states = f(*states)
    new_order = numpy.argsort(f_states)
    basis_permutation = indices[new_order]
    return basis_permutation    

"""
# compute change of basis so that fast_matrix will be block diagonal
# to do this we index states wrt their x[0], x[2], x[1] ordering

np = full_model['np']
state_ordering = lambda *x : (numpy.product(np[1:])*x[0]
                              + numpy.product(np[1])*x[2]
                              + x[1])

block_diag_basis = create_ordered_basis(fast_model, state_ordering)
beta, beta_inverse = create_change_of_basis_matrices(block_diag_basis)

import pylab
pylab.figure()
pylab.spy(beta.todense())
pylab.figure()
pylab.spy(beta_inverse.todense())
pylab.show()
"""