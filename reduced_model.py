import numpy

import rank_approx
import perturb
import cmepy.new_core.cme_solver as cme_solver
import cmepy.core.matrix_cme as matrix_cme
from change_of_basis_snippets import create_ordered_basis

import pylab

def main():
    initial_s_count = 15
    initial_e_count = 6
    epsilon = 0.01
    full_model, slow, fast = rank_approx.get_catalyser_model(initial_s_count,
                                                             initial_e_count,
                                                             epsilon)
    reaction_subsets = (slow, fast)
    (slow_model, fast_model) = perturb.decompose_model(full_model,
                                                       reaction_subsets)
    
    fast_flux_data = cme_solver.create_flux_data(fast_model)
    fast_matrix = matrix_cme.gen_sparse_matrix(full_model['np'],
                                               fast_flux_data)
    fast_matrix *= epsilon
    
    
    reduced_fast_model = {'propensities' : (lambda *x: epsilon*x[1], ),
                          'offset_vectors' : ((1, -1), ),
                          'np' : full_model['np'][1:]}
    
    reduced_fast_flux_data = cme_solver.create_flux_data(reduced_fast_model)
    
    # xxx todo dirty hack of a way to sort by x[0]+x[1]
    # and then break any remaining ties with x[0]-x[1]
    f = lambda *x : (x[0]+x[1])*1000 + (x[0]-x[1])
    states = [numpy.ravel(i) for i in numpy.indices(reduced_fast_model['np'])]
    f_states = f(*states)
    permutation = numpy.argsort(f_states)
    # important! we really need the inverse to plug into the gen_sparse_matrix
    # function
    inverse_permutation = numpy.argsort(permutation)
    
    pylab.figure()
    pylab.plot(inverse_permutation, '.')
    
    reduced_fast_matrix = matrix_cme.gen_sparse_matrix(reduced_fast_model['np'],
                                                       reduced_fast_flux_data,
                                                       inverse_permutation)
    reduced_fast_matrix *= epsilon
    
    pylab.figure()
    pylab.spy(fast_matrix.todense())
    pylab.title('fast matrix')
    pylab.figure()
    pylab.spy(reduced_fast_matrix.todense())
    pylab.title('reduced fast matrix')
    pylab.show()
    pylab.close('all')

if __name__=='__main__':
    main()
