"""
Builds CME matrix for dp/dt, broken into terms for each reaction.
"""

import itertools
import numpy
import scipy.sparse

def optimise_csr_matrix(csr_matrix):
    """
    optimise_csr_matrix(csr_matrix) -> csr_matrix
    
    Performs *in place* operations to optimise the sparse csr matrix data.
    """
    # xxx todo profile performance using permutations / subsets of these
    csr_matrix.sum_duplicates()
    csr_matrix.eliminate_zeros()
    csr_matrix.sort_indices()

def non_neg_states(state_array):
    """
    non_neg_states(state_array) -> bool_array
    
    Returns a boolean array of flags corresponding to those states in
    state_array that have no negative coordinate.
    """
    return numpy.logical_and.reduce(state_array >= 0, axis=0)

def gen_reaction_matrices(model,
                          domain_enum,
                          sink,
                          validity_test):
    """
    gen_reaction_matrices(model, domain_enum, sink, validity_test) ->  generator
    
    Generator yielding the sparse matrices for the dp/dt term of each reaction,
    matching the ordering implied by the ordering of the reaction propensity
    functions and offset vectors in the model.
    
    domain_enum : StateEnum instance enumerating the states in the domain
    sink : boolean flag indicating if the reaction matrices should add
        a 'sink' state used to accumulate probability that flows outside
        of the domain. If sink is set to True, the index of the sink state
        is chosen to be domain_enum.size
    validity_test : validity_test(state_array) -> bool_array
        Returns a boolean array of flags corresponding to those states in
        state_array that are valid.
        
        See: non_neg_states(state_array)
    """
    
    if domain_enum.offset != 0:
        raise NotImplementedError('non-zero domain_enum offset unsupported')
    
    sink = bool(sink)
    if sink:
        sink_index = domain_enum.size
    
    propensities = model['propensities']
    offset_vectors = model['offset_vectors']
    reactions = itertools.izip(propensities, offset_vectors)
    
    src_states = numpy.array(domain_enum.sorted_states)
    src_indices = domain_enum.indices(src_states)
    
    for (propensity, offset_vector) in reactions:
        
        # compute destination states for this offset
        offset_vector = numpy.asarray(offset_vector)[:, numpy.newaxis]
        dst_states = src_states + offset_vector
        
        # determine which states have destination states inside the
        # truncated domain. these will be defined as the 'interior' states.
        # conversely, 'exterior' states are those states of the truncated
        # domain with destination states not in the domain.
        
        interior = domain_enum.contains(dst_states)
        exterior = numpy.logical_not(interior)
        
        num_int_states = numpy.add.reduce(interior)
        num_ext_states = numpy.add.reduce(exterior)
        
        # these lists will be used to accumulate 'COO'-ordinate format
        # sparse matrix data for this reaction.
        
        data = []
        rows = []
        cols = []
        
        # account for the sparse matrix data for the flux out of the
        # interior states of the truncated domain
        
        if num_int_states > 0:
            int_src_indices = numpy.array(src_indices[interior])
            int_src_states = numpy.array(src_states[:, interior])
            int_dst_indices = domain_enum.indices(int_src_states)
            int_coefficients = propensity(*int_src_states)
            
            # flux out
            data.append(-int_coefficients)
            cols.append(int_src_indices)
            rows.append(int_src_indices)
            # flux in
            data.append(int_coefficients)
            cols.append(int_src_indices)
            rows.append(int_dst_indices)
            
        # account for the sparse matrix data for the flux out of the interior
        # states of the truncated domain and into the sink state
        
        if sink and (num_ext_states > 0):
            valid = validity_test(dst_states[:, exterior])
            num_valid_states = numpy.add.reduce(valid)
            if num_valid_states > 0:
                ext_src_indices = numpy.array(src_indices[exterior][valid])
                ext_src_states = numpy.array(src_states[:, exterior][:, valid])
                ext_coefficients = propensity(*ext_src_states)
                
                shape = numpy.shape(ext_src_indices)
                ext_dst_indices = sink_index * numpy.ones(shape,
                                                          dtype = numpy.int)
                
                # these terms account for the flux out of the truncated
                # domain into the sink state
                data.append(-ext_coefficients)
                cols.append(ext_src_indices)
                rows.append(ext_src_indices)
                
                # these terms account for the flux in to the sink state
                # from the truncated domain
                data.append(ext_coefficients)
                cols.append(ext_src_indices)
                rows.append(ext_dst_indices)
        
        matrix_size = domain_enum.size
        if sink:
            matrix_size += 1
        
        matrix_shape = (matrix_size, )*2
        if len(data)==0:
            reaction_matrix = scipy.sparse.csr_matrix(matrix_shape)
        else:            
            # merge data, rows, cols
            data = numpy.concatenate(data)
            cols = numpy.concatenate(cols)
            rows = numpy.concatenate(rows)
            
            # create coo matrix
            coo_data = (data, (rows, cols))
            reaction_matrix = scipy.sparse.coo_matrix(coo_data, matrix_shape)
            
            # convert to sparse csr format, then compress & optimise the storage
            reaction_matrix = reaction_matrix.tocsr()
            optimise_csr_matrix(reaction_matrix)
        
        yield reaction_matrix
    return

def create_diff_eqs(reaction_matrices, phi = None):
    """
    create_diff_eqs(reaction_matrices [, phi]) -> diff_eqs
    
    where diff_eqs(t, p) -> dp_dt
    
    reaction_matrices : sequence of terms of dp/dt matrix corresponding to
        the reactions.
    phi : mapping of time dependent coefficient functions keyed by subsets of
        reaction indices. By default, no time dependent coefficient functions
        are specified, so the returned diff_eqs function is time independent.
    """
    
    num_matrices = len(reaction_matrices)
    if num_matrices == 0:
        raise ValueError('there must be at least one reaction matrix')
    matrix_shapes = set(matrix.shape for matrix in reaction_matrices)
    if len(matrix_shapes) != 1:
        raise ValueError('reaction matrix shapes must all agree')
    matrix_shape = matrix_shapes.pop()
    if len(matrix_shape) != 2:
        raise ValueError('reaction matrices must be two-dimensional')
    if matrix_shape[0] != matrix_shape[1]:
        raise ValueError('reaction matrices must be square')
    
    if phi is None:
        phi = {}
    for reaction_subset in phi:
        if len(reaction_subset) == 0:
            raise ValueError('subsets of reaction indices must be non-empty')
        for i in reaction_subset:
            if not (0 <= i < num_matrices):
                raise ValueError('invalid reaction index: %s' % str(i))
    
    def sum_reaction_matrices(reaction_indices):
        """
        sum_reation_matrices(reaction_indices) -> sum_matrix
        """
        sum_matrix = sum(reaction_matrices[i] for i in reaction_indices)
        optimise_csr_matrix(sum_matrix)
        return sum_matrix
    
    term = {}
    const_indices = set(xrange(num_matrices))
    for reaction_subset in phi:
        const_indices.difference_update(reaction_subset)
        term[reaction_subset] = sum_reaction_matrices(reaction_subset)
    const_indices = frozenset(const_indices)
    if const_indices:
        term[const_indices] = sum_reaction_matrices(const_indices)
    
    def diff_eqs(t, p):
        """
        returns dp / dt for given t and p
        """
        return sum(phi[s](t)*term[s]*p if s in phi else term[s]*p for s in term)
        
    return diff_eqs