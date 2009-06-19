"""
A few utility routines. 27/5/09. R. Fletcher-Costin, ANU
"""

def indices_ext(shape, slices, origin=None):
    """
    An interface to numpy's mgrid routine, supporting simpler slicing notation
    
    Arguments:
    
    shape     : tuple of positive integers giving dimensions of array
    
    slices    : tuple of (abstract) slices indicating which indices we want
                from each dimension
    
    origin    : (optional) tuple of integers giving origin, used to offset
                the returned indices
    """
    
    from itertools import izip
    from numpy import mgrid
    
    assert (len(shape)==len(slices))
    
    # derive concrete slices from abstract slices using provided shape
    slices_concrete = [slice(*sl.indices(n)) for (sl, n) in izip(slices, shape)]
    indices = mgrid[slices_concrete]
    
    if origin is not None:
        assert(len(origin)==len(shape))
        # offset indices by origin
        indices = [x+o for (x, o) in izip(indices, origin)]
    
    return indices 