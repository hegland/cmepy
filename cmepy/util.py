"""
A few utility routines. 27/5/09. R. Fletcher-Costin, ANU
"""

import itertools
import numpy

def consecutive_pairs(p):
    """
    consecutive_pairs(p) -> (p0, p1), (p1, p2), ...
    
    where p is an iterator
    """
    iter = p.__iter__()
    try:
        previous = iter.next()
        while True:
            current = iter.next()
            yield (previous, current)
            previous = current
    except StopIteration:
        return
    
def non_neg(x):
    """
    non_neg(x) -> max(x, 0) [array operation]
    """
    
    return numpy.maximum(x, 0)

def indices_ext(shape, slices=None, origin=None):
    """
    indices_ext(shape [, slices [, origin]]) -> indices array
    
    An interface to numpy's mgrid routine, supporting simpler slicing notation.
    
    Arguments:
    
    shape     : tuple of positive integers giving dimensions of array
    
    slices    : (optional) tuple of (abstract) slices indicating which indices
                we want from each dimension
    
    origin    : (optional) tuple of integers giving origin, used to offset
                the returned indices
    """
    
    if slices is None:
        # include all indices in all dimensions
        slices = (slice(None),)*len(shape)
    
    assert len(shape) == len(slices)
    
    # derive concrete slices from abstract slices using provided shape
    # (mgrid does not accept abstract slices)
    slices_and_dims = itertools.izip(slices, shape)
    slices_concrete = [slice(*sl.indices(n)) for (sl, n) in slices_and_dims]
    indices = numpy.mgrid[slices_concrete]
    
    # only shift indices by origin if origin is non-zero
    # this handles both explicit and implicit (None) zero arguments
    zero_origin = (0, )*len(shape)
    if (origin is not None) and (origin != zero_origin):
        assert len(origin) == len(shape)
        # offset indices by origin
        origin_slice = (slice(None), ) + (numpy.newaxis, )*len(origin)
        indices += numpy.asarray(origin)[origin_slice]
    
    return indices
