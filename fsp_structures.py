import itertools
import numpy
import scipy.sparse
import lexarrayset

class StateEnumeration(object):
    """
    Maintains a state, index bijection
    """
    def __init__(self, dim):
        """
        Initialise the state enumeration using the provided
        LexArraySet initial_states.
        """
        self.dim = dim
        self.ordered_states = None
        self.unordered_states = None
        self.index = None
        self.size = None
        self.offset = 0
        self.reinitialise(lexarrayset.empty(self.dim))
    
    def update_ordering(self):
        """
        Used to maintain the ordering invariants.
        """
        order = numpy.lexsort(self.unordered_states)
        self.index = self.index[order]
        self.size = numpy.size(self.index)
        if self.size > 0:
            self.ordered_states = self.unordered_states[:, order]
        else:
            self.ordered_states = self.unordered_states
    
    def extend(self, sigma):
        """
        Adds the states in the LexArraySet sigma to the state enumeration.
        
        These states must be disjoint to the existing states in this
        enumeration.
        
        The indexing of existing states in this enumeration will be
        unchanged.
        """
        
        self.unordered_states = numpy.hstack((self.unordered_states,
                                              sigma.data))
        index = numpy.arange(self.size, self.size+sigma.data.shape[1])
        self.index = numpy.concatenate((self.index, index))
        self.update_ordering()
    
    def reinitialise(self, initial_states):
        """
        reinitialise the StateEnumeration with the given initial states
        """
        self.unordered_states = initial_states.data
        self.index = numpy.arange(self.unordered_states.shape[1])
        self.update_ordering()
        self.offset = 0
        
    def indices(self, states):
        """
        returns an array of the enumeration indices for the
        states stored in the LexArraySet states.
        
        If states is not a LexArraySet, ie if states is a 2d array of 
        possible non-unique rows, attempt to return the right thing anyhow ...
        """
        if type(states) is lexarrayset.LexArraySet:
            index = numpy.array(self.index[lexarrayset.member(self.ordered_states,
                                                              states.data)])
            return index + self.offset
        else:
            # assume states is a two dimensional array with
            # potentially non unique rows
            
            # due to non-uniqueness, there is a bit of messing
            # about in order to reduce the states to
            # a unique set, find the indices for those states, then
            # invert the unique operation
            
            print 'StateEnumeration.indices : non-unique case ...'
            
            unique_states, inverse = lexarrayset.unique(states,
                                                        return_inverse=True)
            las = lexarrayset.create(unique_states)
            las_indices = self.indices(las)
            return las_indices[inverse]
    
    def states(self, index):
        """
        returns an array of the states corresponding to the
        provided enumeration index array.
        """
        
        return numpy.array(self.unordered_states[:, index - self.offset])
    
    def pack_distribution(self, p, y=None):
        """
        convenience routine to translate a distribution from a dictionary to
        a dense array, using this state enumeration 
        """
        
        if y is None:
            y = numpy.zeros((self.size, ), dtype=numpy.float)
        
        states, values = zip(*p.items())
        # convert from list of states to list of coordinate vectors
        states = numpy.array(zip(*states))
        
        # guard against case where p is empty
        if len(states) == 0:
            return y
        # now sort the states, keeping them synchronised with the
        # ordering of the values
        order = numpy.lexsort(states)
        states = lexarrayset.create(states[:, order], unique_data = True)
        values = numpy.array(values)[order]
        indices = self.indices(states)
        y[indices] = values
        return y
        
    
    def unpack_distribution(self, y, p=None):
        """
        convenience routine to translate a distribution from a dense array
        to a dictionary, using this state enumeration
        """
        indices = numpy.arange(numpy.size(y))
        states = self.states(indices)
        # convert from list of coordinate vectors to list of states
        states = zip(*states)
        if p is None:
            p = {}
        for state, value in itertools.izip(indices, states):
            if value != 0.0:
                p[state] = value
        return p
    
class ShiftStateUpdate(object):
    def __init__(self,
                 delta_interior,
                 delta_minus_boundary,
                 delta_plus_boundary):
        self.delta_interior = delta_interior
        self.delta_minus_boundary = delta_minus_boundary
        self.delta_plus_boundar = delta_plus_boundary

class ShiftStateTracker(object):
    """
    this thing tracks different subsets of source states
    for a specific state space shift
    """
    
    def __init__(self, shift, dim):
        self.shift = shift
        self.dim = dim
        self.interior = lexarrayset.empty(self.dim)
        self.boundary = lexarrayset.empty(self.dim)
    
    def extend(self, sigma):
        """
        sigma should be a lexical array set of the new states
        for the state space
        """
        
        shift_sigma = sigma.shift(self.shift)
        shift_sigma_inverse = sigma.shift(-self.shift)
        
        # nb scheme here supposes that sigma is fairly small compared
        # to interior and boundary, and so only applies shifts
        # to sigma (or derived values), not to interior or boundary
        
        # find states in boundary with destinations in sigma
        bdry_to_sig, bdry_to_sink = self.boundary.split(shift_sigma_inverse)
        # find states in sigma with destinations in interior
        sig_to_int, sig_to_int_c = shift_sigma.split(self.interior)
        sig_to_int.shift_update(-self.shift)
        # find states in sigma with destinations in boundary
        sig_to_bdry, sig_to_bdry_c = sig_to_int_c.split(self.boundary)
        sig_to_bdry.shift_update(-self.shift)
        # find states in sigma with destinations in sigma
        sig_to_sig, sig_to_sink = sig_to_bdry_c.split(sigma)
        sig_to_sig.shift_update(-self.shift)
        sig_to_sink.shift_update(-self.shift)
        
        delta_interior = reduce(lambda x, y : x.union(y),
                                (bdry_to_sig,
                                 sig_to_int,
                                 sig_to_bdry,
                                 sig_to_sig))
        
        delta_minus_boundary = self.boundary.difference(bdry_to_sink)
        delta_plus_boundary = sig_to_sink
        
        self.interior.union_update(delta_interior)
        self.boundary.difference_update(delta_minus_boundary)
        self.boundary.union_update(delta_plus_boundary)
        
        assert self.interior.intersection(self.boundary).size == 0
        
        update = ShiftStateUpdate(delta_interior,
                                  delta_minus_boundary,
                                  delta_plus_boundary)
        return update

def create_flux_matrix(source_state_enum,
                       source_states,
                       shift,
                       propensity,
                       matrix_size,
                       projections=None):
    
    source_indices = source_state_enum.indices(source_states)
    propensities = propensity(source_states)
    
    row = []
    col = []
    val = []
    
    row.append(source_indices)
    col.append(source_indices)
    val.append(-propensities)
    
    dest_states = source_states.shift(shift)
    if projections is None:
        dest_indices = source_state_enum.indices(dest_states)
        row.append(dest_indices)
        col.append(source_indices)
        val.append(propensities)
    else:
        for (dest_state_enum, transform) in projections:
            dest_states = transform(dest_states.data)
            dest_indices = dest_state_enum.indices(dest_states)
            row.append(dest_indices)
            col.append(source_indices)
            val.append(propensities)
    
    row = numpy.concatenate(row)
    col = numpy.concatenate(col)
    val = numpy.concatenate(val)
        
    coo_data = ((row, col), val)
    coo_shape = (matrix_size )*2
    coo = scipy.sparse.coo_matrix(coo_data, coo_shape)
    csr = coo.tocsr()
    # xxx todo profile this, see if it actually helps
    csr.sum_duplicates()
    csr.eliminate_zeros()
    csr.sort_indices()
    return csr

class FluxMatrixSum(object):
    def __init__(self):
        self.sum = None
    
    def add_term(self, term):
        if self.sum is None:
            self.sum = term
        else:
            self.sum += term


class InteriorMatrixFactory(object):
    """
    Wrapper for create_flux_matrix routine for producing
    interior flux matrices
    """
    def __init__(self, state_enum, shift, propensity):
        self.state_enum = state_enum
        self.shift = shift
        self.propensity = propensity
    
    def produce_term(self, source_states, matrix_size):
        term = create_flux_matrix(self.state_enum,
                                  self.state_enum,
                                  source_states,
                                  self.shift,
                                  self.propensity,
                                  matrix_size)
        return term

class BoundaryMatrixFactory(object):
    """
    Wrapper for create_flux_matrix routine for producing
    boundary flux matrices
    """
    
    def __init__(self,
                 source_state_enum,
                 dest_state_enum,
                 shift,
                 propensity,
                 dest_state_transform):
        self.source_state_enum = source_state_enum
        self.dest_state_enum = dest_state_enum
        self.shift = shift
        self.propensity = propensity
        self.dest_state_transform = dest_state_transform
    
    def produce_term(self, source_states, matrix_size):
        term = create_flux_matrix(self.source_state_enum,
                                  self.dest_state_enum,
                                  source_states,
                                  self.shift,
                                  self.propensity,
                                  matrix_size,
                                  self.dest_state_transform)
        return term

class DiffEqsFactory(object):
    
    def __init__(self, reactions_to_time_deps, time_deps):
        self.reactions_to_time_deps = reactions_to_time_deps
        self.time_deps = time_deps
        self.interior_matrices = {}
        self.boundary_matrices = {}
    
    
    def update_boundary_matrix(self, reaction, boundary):
        self.boundary_matrices[reaction] = boundary
        
    def update_interior_matrix(self, reaction, interior):
        time_dep = self.reactions_to_time_deps[reaction]
        if time_dep not in self.interior_matrices:
            self.interior_matrices[time_dep] = interior
        else:
            current_interior = self.interior_matrices[time_dep]
            if current_interior.shape == interior.shape:
                self.interior_matrices[time_dep] += interior
            else:
                shape = tuple(numpy.maximum(current_interior.shape,
                                            interior.shape))
                
                if shape != current_interior.shape:
                    current_interior.set_shape(shape)
                if shape != interior.shape:
                    interior.set_shape(shape)
                self.interior_matrices[time_dep] = current_interior + interior
    
        
    
    def create_diff_eqs(self):
        self.matrices = {}
        for time_dep in self.interior_matrices:
            interior = self.interior_matrices[time_dep]
            self.matrices[time_dep] = interior
        for reaction in self.boundary_matrices:
            time_dep = self.reactions_to_time_deps[reaction]
            boundary = self.boundary_matrices[reaction]
            # nb dont use += as we do not wish to alter the matrix
            # stored by reference in matrices
            self.matrices[time_dep] = self.matrices[time_dep] + boundary
        
        def diff_eqs(t, p):
            dp_dt = None
            for time_dep in self.matrices:
                term = self.matrices[time_dep]*p
                phi = self.time_deps[time_dep]
                if phi is not None:
                    term *= phi(t)
                if dp_dt is None:
                    dp_dt = term
                else:
                    dp_dt += term
            return dp_dt
        
        return diff_eqs
        
        