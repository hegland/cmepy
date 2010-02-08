import itertools

import numpy

import cmepy.new_core.domain as domain
import cmepy.new_core.lexarrayset as lexarrayset


def create(initial_states):
    """
    create(initial_states) -> StateEnum instance
    
    instantiates a StateEnum instance using the provided 'initial_states'.
    """
    return StateEnum(initial_states)

class StateEnum(object):
    """
    Maintains bijection between set of n unique states and range(n)
    """
    def __init__(self, initial_states):
        """
        Initialise the state enumeration using the provided
        'initial_states'.
        """
        self.ordered_states = None
        self.unordered_states = None
        self.index = None
        self.size = None
        self.offset = 0
        self.reinitialise(initial_states)
    
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
        Adds the states in the array 'sigma' to the state enumeration.
        
        These states must be disjoint to the existing states in this
        enumeration.
        
        The indexing of existing states in this enumeration will be
        unchanged.
        """
        
        sigma_unique = lexarrayset.unique(numpy.asarray(sigma))
        
        self.unordered_states = numpy.hstack((self.unordered_states,
                                              sigma_unique))
        index = numpy.arange(self.size, self.size+sigma_unique.shape[1])
        self.index = numpy.concatenate((self.index, index))
        self.update_ordering()
    
    def reinitialise(self, initial_states):
        """
        reinitialise the StateEnumeration with the given 'initial_states'
        """
        
        initial_states = numpy.asarray(initial_states)
        
        self.unordered_states = lexarrayset.unique(initial_states)
        self.index = numpy.arange(self.unordered_states.shape[1])
        self.update_ordering()
        self.offset = 0
        
    def indices(self, states):
        """
        indices(states) -> index_array
        
        returns an array of the enumeration indices for the
        states stored in the array 'states'.
        """
        
        states = numpy.asarray(states)
        
        # assume states is a two dimensional array with
        # potentially non unique rows
        
        # due to non-uniqueness, there is a bit of messing
        # about in order to reduce the states to
        # a unique set, find the indices for those states, then
        # invert the unique operation
        
        unique_states, unique_inverse = lexarrayset.unique(states,
                                                           return_inverse=True)
        
        # subtlety : we need the boolean array members to correspond
        # to the ordered states and thus also to the current index,
        # hence we test to see which elements of ordered_states
        # are contained in the unique states
        #
        # note that this differs from the members query in the
        # contains method
        
        members = lexarrayset.member(self.ordered_states, unique_states)
        member_index = numpy.array(self.index[members] + self.offset)
        return member_index[unique_inverse]
    
    def contains(self, states):
        """
        contains(states) -> bool_array
        
        returns a boolean array of flags indicates which of the
        states stored in the array 'states' are contained in the
        state enumeration.
        """
        
        states = numpy.asarray(states)
        
        unique_states, unique_inverse = lexarrayset.unique(states,
                                                           return_inverse=True)
        
        # subtlety : we need the boolean array members to correspond to the
        # unique states, hence we test to see which elements of the unique
        # states are contained in the ordered states
        #
        # note that this differs from the members query in the indices
        # method
        
        members = lexarrayset.member(unique_states, self.ordered_states)
        return members[unique_inverse]
        
    
    def states(self, index):
        """
        returns an array of the states corresponding to the
        provided enumeration index array.
        """
        
        index = numpy.asarray(index)
        
        return numpy.array(self.unordered_states[:, index - self.offset])
    
    def pack_distribution(self, p_sparse, p_dense=None):
        """
        convenience routine to translate a distribution from a dictionary to
        a dense array, using this state enumeration 
        """
        
        if p_dense is None:
            p_dense = numpy.zeros((self.size, ), dtype=numpy.float)
        
        # guard against case where p_sparse is empty
        if len(p_sparse) == 0:
            return p_dense
        
        p_states, p_values = domain.from_mapping(p_sparse)
        
        # now sort the states, keeping them synchronised with the
        # ordering of the values
        order = numpy.lexsort(p_states)
        p_states = p_states[:, order]
        p_values = p_values[order]
        p_indices = self.indices(p_states)
        p_dense[p_indices] = p_values
        return p_dense
        
    
    def unpack_distribution(self, p_dense, p_sparse=None):
        """
        convenience routine to translate a distribution from a dense array
        to a dictionary, using this state enumeration
        """
        p_indices = numpy.arange(numpy.size(p_dense))
        # convert from list of coordinate vectors to list of states
        p_states = domain.to_iter(self.states(p_indices))
        if p_sparse is None:
            p_sparse = {}
        for index, state in itertools.izip(p_indices, p_states):
            value = p_dense[index]
            if value != 0.0:
                p_sparse[state] = value
        return p_sparse
