"""
Generates sparse matrix representations of CME models
"""

import numpy
import scipy.sparse
from cmepy.core.core_cme import process_offset_vectors, gen_reaction_actions

def gen_sparse_matrix(np, reaction_actions):
    """
    Generates a scipy.sparse.coo_matrix for the given model np and
    reaction_actions.
    """
    
    num_states = numpy.product(np)
    flat_states = numpy.arange(num_states)
    inflated_states = numpy.reshape(flat_states, np)
    
    matrix_rows = numpy.array(())
    matrix_columns = numpy.array(())
    matrix_coefficients = numpy.array(())
    
    for action in reaction_actions:
        coefficients, source_indices, dest_indices = action
        
        # account for flux out of source states
        matrix_columns = numpy.append(matrix_columns,
                                       inflated_states[source_indices])
        matrix_rows = numpy.append(matrix_rows,
                                   inflated_states[source_indices])
        matrix_coefficients = numpy.append(matrix_coefficients,
                                           -1.0*coefficients)
        
        # account for flux in to destination states
        matrix_columns = numpy.append(matrix_columns,
                                       inflated_states[source_indices])
        matrix_rows = numpy.append(matrix_rows,
                                   inflated_states[dest_indices])
        matrix_coefficients = numpy.append(matrix_coefficients,
                                           coefficients)
    
    matrix_data = (matrix_coefficients, (matrix_rows, matrix_columns))
    matrix_shape = (num_states, num_states)
    sparse_matrix = scipy.sparse.coo_matrix(matrix_data, shape = matrix_shape)
    return sparse_matrix
    
    
def sparse_matrix_factory(**args):
    """
    Creates and returns the sparse matrix A representation for the CME
    
    dp/dt = Ap(t) .
    
    (NB since A is a constant, this method does not support models featuring
    time dependent propensity functions)
    """
    
    np = args['np']
    propensities = args['propensities']
    norigin = args.get('norigin', None)
    offset_vectors = args.get('offset_vectors', None)
    offset_vectors = process_offset_vectors(np,
                                            propensities,
                                            offset_vectors) 
    reaction_actions = gen_reaction_actions(np,
                                            propensities,
                                            offset_vectors,
                                            norigin)
    time_dependence = args.get('time_dependence', None)
    if time_dependence is not None:
        raise ValueError('time dependence is not supported by matrices')
    
    return gen_sparse_matrix(np, reaction_actions)

def spare_matrix_diff_eqs_factory(sparse_a_matrix):
    """
    returns the function diff_eqs(t, p) := sparse_a_matrix*p
    """
    
    def diff_eqs(t, p):
        """
        this routine defines our mapping diff_eqs(t,p) ---> dp/dt(t,p)
        which is later called by the ode solver
        """
        return sparse_a_matrix*p
    
    return diff_eqs
        

def cme(sparse_matrix, **args):
    """
    Creates a scipy.integrate.ode instance for the specified cme problem
    
    The underlying ode solver employed is vode, in mode bdf.    
    """
    import numpy
    from scipy.integrate import ode
    np = args['np']
    num_states = numpy.prod(np)
    p_0 = args.get('p0', None)
    t_0 = args.get('t0', 0.0)
    
    if p_0 is None:
        p_0 = numpy.zeros(num_states)
        p_0[0] = 1.0
    else:
        p_0 = numpy.reshape(p_0, (num_states,))
    
    # create the differential equations function
    diff_eq = spare_matrix_diff_eqs_factory(sparse_matrix)
    # create the ode and set the integrator type and method
    cme_de = ode(diff_eq).set_integrator('vode', method='bdf')
    # set the initial values
    cme_de.set_initial_value(p_0, t_0)
    
    return cme_de