"""
Various utility routines
"""

import itertools
import numpy

def shape_invariant(f):
    """
    Returns wrapped copy of f that ensures output shape matches input shape.
    
    Details : if the shapes of the first argument and output of the function f 
    do not agree, the output will be expanded via ``numpy.tile`` until it is
    the same shape as the first input. Shapes are computed via ``numpy.shape``.
    """
    
    def shape_invariant_f(*args):
        x = f(*args)
        in_shape = numpy.shape(args[0])
        out_shape = numpy.shape(x)
        if in_shape != out_shape:
            out_shape += (1, )*(len(in_shape) - len(out_shape))
            tiling = tuple(x / y for (x, y) in zip(in_shape, out_shape))
            x = numpy.tile(x, tiling)
        return  x
    return shape_invariant_f
    
def consecutive_pairs(p):
    """
    Returns (p0, p1), (p1, p2), ...
    
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
    Returns max(x, 0) [array operation]
    """
    
    return numpy.maximum(x, 0)

def indices_ext(shape, slices=None, origin=None):
    """
    Returns an array of indices.
    
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
