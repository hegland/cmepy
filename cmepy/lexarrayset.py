"""
lexical array set operations

these operations are based upon the one dimensional array set operations
from numpy.lib.arraysetops, but generalised to work for sets of m-tuples,
where each element is stored as a row of a 2d m by n array, using numpy's
'lexsort' lexical sorting function.
"""

import numpy

def unique(las, return_inverse=False):
    """
    returns a sorted vector of unique states
    
    if the optional flag return_inverse is set to True,
    additionally returns an index vector used to
    inverse the unique operation and recover the
    original vector
    """
    
    # argsort the array via lexical sorting using the keys
    # las[0, :] to las[-1, :], in increasing priority
    order = numpy.lexsort(las)
    if numpy.size(order) == 0:
        return las
    slas = las[:, order]
    # then figure out the indices of the first instance of each row
    not_equal_adj = numpy.logical_or.reduce(slas[:, :-1] != slas[:, 1:])
    not_equal_adj = numpy.concatenate(([True], not_equal_adj))
    
    uslas = slas[:, not_equal_adj]
    if return_inverse:
        order_inverse = order.argsort(kind='mergesort')
        # compute the unique inverse indices by summing over the
        # boolean array. because we want to compute 0-based indices
        # it is necessary to set the first boolean to False.
        # (alternatively, we could have subtracted 1 from the entire
        # result after summing)
        not_equal_adj[0] = False
        unique_inverse = numpy.add.accumulate(not_equal_adj)[order_inverse]
        return uslas, unique_inverse
    else:
        return uslas

def nonunique_member(arr1, las2):
    """
    vectorised set membership operation for lexical array arr1 and
    lexical array set las2
    
    in general, the rows of array arr1 can be non-unique
    
    returns a boolean array 'mask' such that
    arr1[:, mask] is the subset of rows of arr1 that are also
    rows of las2
    """
    las1, unique_inverse = unique(arr1, return_inverse=True)
    return member(las1, las2)[unique_inverse]

def member(las1, las2):
    """
    vectorised set membership operation for lexical array sets las1, las2
    
    returns a boolean array 'mask' such that
    
    las1[:, mask] is the subset of rows of las1 that
    are also rows of las2
    """
    
    las = numpy.hstack((las1, las2))
    
    las1_n = numpy.shape(las1)[1]
    
    # since the 'in' operation is non-commutative we must
    # use a stable sort algorithm. this ensures that
    # if slas[i] == slas[i+1], then slas[i] is the element
    # from las1, while slas[i+1] is the element from las2
    
    # by grepping through the numpy source it seems that lexsort
    # should be a stable sort (is this officially documented anywhere?)
    
    order = numpy.lexsort(las, )
    
    if numpy.size(order) == 0:
        return numpy.zeros((las1_n, ), dtype=numpy.bool)
    slas = las[:, order]
    equal_adj = numpy.logical_and.reduce(slas[:, :-1] == slas[:, 1:])
    mask = numpy.concatenate((equal_adj, [False]))
    
    inverse_order = order.argsort(kind='mergesort')[:las1_n]    
    return mask[inverse_order]

def split(las1, las2):
    """
    returns (las1 intersect las2, las1 difference las2)
    """
    if numpy.size(las1) == 0:
        # if las1 is empty, return a couple of copies of las1
        return (numpy.array(las1), numpy.array(las1))
    mask = member(las1, las2)
    return (numpy.array(las1[:, mask]),
            numpy.array(las1[:, numpy.logical_not(mask)]))

def difference(las1, las2):
    """
    returns las1 difference las2
    """
    if numpy.size(las1) == 0:
        # if las1 is empty, return a copy of las1
        return numpy.array(las1)
    return numpy.array(las1[:, numpy.logical_not(member(las1, las2))])

def intersection(las1, las2):
    """
    intersection of las1 with las2
    """
    las = numpy.hstack((las1, las2))
    order = numpy.lexsort(las)
    if numpy.size(order) == 0:
        return las
    slas = las[:, order]
    equal_adj = numpy.logical_and.reduce(slas[:, :-1] == slas[:, 1:])
    return slas[:, equal_adj]

def union(las1, las2):
    """
    union of las1 with las2
    """
    return unique(numpy.hstack((las1, las2)))

def empty(dim):
    """
    returns an empty LexArraySet of dimension dim.
    """
    empty_data = numpy.zeros((dim, 0), dtype=numpy.int)
    return LexArraySet(empty_data)

def create(data, unique_data=False):
    """
    returns a new LexArraySet for the given data
    """
    return LexArraySet(data, unique_data)

class LexArraySet(object):
    """
    LexArraySet is an implementation of a set as a 2d array, where
    the members of the set are the rows of the array. The rows
    are ordered using lexical ordering.
    """
    
    def __init__(self, data, unique_data=False):
        """
        data can either be another LexArraySet instance, in which
        case a copy is made of that instance's data, or a
        two-dimensional numpy array, where each row is interpreted
        as a tuple belonging to the set.
        
        If data is a two-dimensional numpy array, then the optional
        unique_data flag may be set to True to indicate that the
        rows of data are already unique.
        """
        if type(data) is LexArraySet:
            self.data = numpy.array(data.data)
        else:
            data = numpy.asarray(data)
            if not unique_data:
                self.data = unique(data)
            else:
                self.data = data
    
    def _get_size(self):
        """
        number of elements in set (equal to number of rows of the lexical array)
        """
        shape = numpy.shape(self.data)
        if len(shape) < 2:
            return 0
        else:
            return shape[1]
    
    size = property(_get_size)
    
    def member(self, rhs):
        """
        las1.member(las2) -> mask; mask[i] True iff row i of las1 is in las2
        """
        return member(self.data, rhs.data)
    
    def split(self, rhs):
        """
        las1.split(las2) -> (las1.intersect(las2), las1.difference(las2))
        """
        intersect, diff = split(self.data, rhs.data)
        las_intersect = LexArraySet(intersect, unique_data=True)
        las_diff = LexArraySet(diff, unique_data=True)
        return las_intersect, las_diff
    
    def difference(self, rhs):
        """
        las1.difference(las2) -> diff; diff's rows are those of las1 not in las2
        """
        return LexArraySet(difference(self.data, rhs.data), unique_data=True)
    
    def intersection(self, rhs):
        """
        las1.intersection(las2) -> isect; isect's rows common to las1, las2
        """
        return LexArraySet(intersection(self.data, rhs.data), unique_data=True)
    
    def union(self, rhs):
        """
        las1.union(las2) -> u; u's rows are union of rows in las1, las2
        """
        return LexArraySet(union(self.data, rhs.data), unique_data=True)
    
    def shift(self, offset):
        """
        las.shift(offset) -> slas; where rows of slas are rows of las + offset
        
        offset must be of compatible shape to the rows of las.
        """
        offset = numpy.asarray(offset)[:, numpy.newaxis]
        return LexArraySet(self.data + offset, unique_data=True)
    
    def difference_update(self, rhs):
        """
        in place difference
        """
        self.data = difference(self.data, rhs.data)
    
    def intersection_update(self, rhs):
        """
        in place intersection
        """
        self.data = intersection(self.data, rhs.data)
    
    def union_update(self, rhs):
        """
        in place union
        """
        self.data = union(self.data, rhs.data)
    
    def shift_update(self, shift):
        """
        in place shift
        """
        shift = numpy.asarray(shift)[:, numpy.newaxis]
        self.data += shift
