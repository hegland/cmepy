import itertools
import numpy
import scipy.sparse

class StateIndexMap(object):
    def __init__(self, dim, vectstates, origin=None):
        self.dim = dim
        self.vectstates = numpy.asarray(vectstates)
        if origin is None:
            self.origin = 0
        else:
            self.origin = int(origin)
        
        # build dict-based state<->index transforms
        self.state_to_index = {}
        self.index_to_state = {}
        for relative_index, state in enumerate(itertools.izip(*vectstates)):
            index = self.origin + relative_index
            self.state_to_index[state] = index
            self.index_to_state[index] = state
        
        self.size = len(self.state_to_index)
        
        # dynamically generate numpy-vectorised state<->index transforms
        v_func = numpy.vectorize(lambda *state : self.state_to_index[state],
                                 [numpy.int])
        self.vectstates_to_indices = lambda states : v_func(*states)
        
        func = lambda index : self.index_to_state[index]
        otypes = (numpy.int, )*self.dim   
        self.indices_to_vectstates = numpy.vectorize(func, otypes)
        
        # generate characteristic function for states
        func2 = lambda *state : state in self.state_to_index
        v_func2 = numpy.vectorize(func2, [numpy.bool])
        self.vectstates_in_map = lambda states : v_func2(*states)
        
        # create constants
        self.indices = self.vectstates_to_indices(self.vectstates)
    
    def pack_distribution(self, p, y=None):
        if y is None:
            y = numpy.zeros((self.size, ))
        for state in p:
            index = self.state_to_index[state] - self.origin
            y[index] = p[state]
        return y
    
    def unpack_distribution(self, y, p=None):
        if p is None:
            p = {}
        for index in self.index_to_state:
            mass = y[index - self.origin]
            if mass != 0.0:
                state = self.index_to_state[index]
                p[state] = mass
        return p
            
        
def create_flux_matrices(model,
                         state_index_map,
                         error_trackers,
                         error_projection_dim):
    
    propensities = model['propensities']
    offset_vectors = model['offset_vectors']
    
    source_states = state_index_map.vectstates
    source_indices = state_index_map.indices
    
    flux_matrices = []
    for (propensity, offset_vector) in itertools.izip(propensities,
                                                      offset_vectors):
        
        # ensure offset_vector has the correct structure
        offset_vector = numpy.asarray(offset_vector)[:, numpy.newaxis]
        dest_states = source_states + offset_vector
        
        # figure out which of these states are legal
        state_mask = state_index_map.vectstates_in_map(dest_states)
        not_state_mask = numpy.logical_not(state_mask)
        
        num_states = numpy.add.reduce(state_mask)
        num_error_states = numpy.add.reduce(not_state_mask)
        
        data = []
        rows = []
        cols = []
        
        if num_states > 0:
            masked_source_indices = numpy.array(source_indices[state_mask])
            masked_source_states = numpy.array(source_states[:, state_mask])
            masked_dest_indices = state_index_map.vectstates_to_indices(dest_states[:, state_mask])
            coefficients = propensity(*masked_source_states)
            
            data.append(-coefficients)
            cols.append(masked_source_indices)
            rows.append(masked_source_indices)
            data.append(coefficients)
            cols.append(masked_source_indices)
            rows.append(masked_dest_indices)
            
        
        if num_error_states > 0:
            # only operate on valid states that are not contained in the truncated state space
            validity_mask = numpy.logical_and.reduce(dest_states[:, not_state_mask]>=0, axis=0)
            num_valid_states = numpy.add.reduce(validity_mask)
            
            if num_valid_states > 0:
                valid_dest_states = numpy.array(dest_states[:, not_state_mask][:, validity_mask])
                valid_source_indices = numpy.array(source_indices[not_state_mask][validity_mask])
                valid_source_states = numpy.array(source_states[:, not_state_mask][:, validity_mask])
                coefficients = propensity(*valid_source_states)
                
                # flux lost due to truncation of state space
                data.append(-coefficients)
                cols.append(valid_source_indices)
                rows.append(valid_source_indices)
                
                # record lost flux in error states, if any are provided
                for error_tracker in error_trackers:
                    error_projection, error_index_map = error_tracker
                    
                    # xxx todo this doesn't seem to work properly if the
                    # error projection dim is set to 1
                    v = numpy.vectorize(error_projection, otypes = [numpy.int, ]*error_projection_dim)
                    
                    dest_error_states = v(*valid_dest_states)
                    
                    dest_error_indices = error_index_map.vectstates_to_indices(dest_error_states)
                    
                    data.append(coefficients)
                    cols.append(valid_source_indices)
                    rows.append(dest_error_indices)
        
        # figure out how many indices there are in total
        size = state_index_map.size
        for error_tracker in error_trackers:
            size += error_tracker[1].size
        
        if len(data)==0:
            flux_matrix = scipy.sparse.csr_matrix((size, )*2)
        else:
            def join(x):
                if len(x) == 0:
                    return []
                
                size = numpy.add.reduce([numpy.size(a) for a in x])
                joined_x = numpy.zeros((size, ), dtype=x[0].dtype)
                offset = 0
                for a in x:
                    a_size = numpy.size(a)
                    joined_x[offset:offset+a_size] = a
                    offset += a_size
                return joined_x
            
            # merge data, rows, cols
            data = join(data)
            cols = join(cols)
            rows = join(rows)
            
            flux_matrix = scipy.sparse.coo_matrix((data, (rows, cols)), (size, )*2)
            #print '\t coo : %s' % repr(flux_matrix)
            flux_matrix = flux_matrix.tocsr()
            #print '\t csr : %s' % repr(flux_matrix)
        flux_matrix.sum_duplicates()
        flux_matrix.eliminate_zeros()
        #print '\t csr (cleaned?) : %s' % repr(flux_matrix)
        flux_matrices.append(flux_matrix)
            
    return flux_matrices

def create_diff_eqs(size, flux_matrices, time_dependencies = None):
    """
    returns diff_eqs(t, y) instance to pass to Solver
    
    size : length of y vector
    flux_matrices : list of flux_matrices for reactions
    time_dependencies : dict of time dependent coefficient functions
        keyed by subsets of reaction indices
    """
    
    if time_dependencies is None:
        time_dependencies = {}
    assert type(time_dependencies) is dict
    
    def sum_flux_matrices(reaction_indices):
        sum_matrix = scipy.sparse.csr_matrix((size, )*2)
        for reaction_index in reaction_indices:
            sum_matrix = sum_matrix + flux_matrices[reaction_index]
        sum_matrix.sum_duplicates()
        sum_matrix.eliminate_zeros()
        return sum_matrix
    
    matrices = {}
    const_indices = set(range(len(flux_matrices)))
    for reaction_subset in time_dependencies:
        const_indices.difference_update(reaction_subset)
        matrices[reaction_subset] = sum_flux_matrices(reaction_subset)
    
    const_indices = frozenset(const_indices)
    if const_indices:
        matrices[const_indices] = sum_flux_matrices(const_indices)
    
    assert len(matrices) > 0, 'there must be at least one reaction'
    
    def diff_eqs(t, p):
        """
        returns dp / dt for given t and p
        """
        dp_dt = None
        for reaction_subset in matrices:
            term = matrices[reaction_subset]*p
            if reaction_subset in time_dependencies:
                term = term * time_dependencies[reaction_subset](t)
            if dp_dt is not None:
                dp_dt += term
            else:
                dp_dt = term
        return dp_dt
        
    return diff_eqs