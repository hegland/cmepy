"""
test aggregation routines for cme singlular perturbation calculations
"""

import numpy
import scipy.sparse

def extract_sparse_blocks(a):
    """
    determines the block-diagonal structure of a
    
    a should be a sparse matrix
    
    returns a list of entries of the form
    
    (block_start, block_size, block_data)
    
    where block_start gives the index of the first element
    of a along the diagonal contained by the block, block_size gives
    the length of the block (along the diagonal) and block_data
    is a sparse matrix containing the block entries themselves
    (with (row, col) coordinates translated by (-block_size, -block_size))
    """
    
    # translate to sparse COOrdinate format
    coo_a = a.tocoo()
    
    row = numpy.array(coo_a.row)
    col = numpy.array(coo_a.col)
    data = numpy.array(coo_a.data)
    
    # remove zero entries
    nonzero_indices = numpy.nonzero(data)[0]
    row = row[nonzero_indices]
    col = col[nonzero_indices]
    data = data[nonzero_indices]
    
    num_entries = numpy.size(data)
    
    if num_entries == 0:
        return []
    
    # min and max coords
    min_coord = numpy.minimum(row, col)
    max_coord = numpy.maximum(row, col)
    
    # argsort all entry data with respect to increasing min coord
    min_coord_order = numpy.argsort(min_coord)
    
    row = row[min_coord_order]
    col = col[min_coord_order]
    data = data[min_coord_order]
    min_coord = min_coord[min_coord_order]
    max_coord = max_coord[min_coord_order]
    
    # process entries one by one
    # idea : use a moving window along the diagonal,
    # from block_start to block_end. entries contribute to the block
    # containing the interval [min_coord, max_coord] along the diagonal
    
    block_start = 0
    block_end = 0
    j = 0
    blocks = []

    for i in xrange(num_entries):
        entry_start = min_coord[i]
        entry_end = max_coord[i]
        
        if entry_start <= block_end:
            if entry_end > block_end:
                block_end = entry_end
        else:
            if i > j:
                # finalise previous block
                block_data = (data[j:i],
                              (row[j:i]-block_start,
                               col[j:i]-block_start))
                blocks.append((block_start,
                               block_end-block_start+1,
                               scipy.sparse.coo_matrix(block_data)))
            # commence new block
            block_start = entry_start
            block_end = entry_end
            j = i
    
    # finalise last block
    block_data = (data[j:],
                  (row[j:]-block_start,
                   col[j:]-block_start))
    blocks.append((block_start,
                   block_end-block_start+1,
                   scipy.sparse.coo_matrix(block_data)))
    
    return blocks

def compute_aggregation(a):
    blocks = extract_sparse_blocks(a)


def print_sparse_blocks(blocks):
    for i, (start, size, data) in enumerate(blocks):
        print '** block %d' % i
        print '\tstart %d' % start
        print '\tsize %d' % size
        print '\tdata:'
        print str(data)
        print ''

def display_sparse_blocks(blocks):
    import pylab
    import itertools
    pylab.figure()
    for start, size, data in blocks:
        pylab.scatter(data.col+start,
                      data.row+start)
    pylab.show()

def test_extract_sparse_blocks_a():
    a_dense = numpy.zeros((10, 10))
    a_dense[0, 0] = 1
    a_dense[1, 1] = 1
    a_dense[1, 2] = 1
    a_dense[2, 2] = 1
    a = scipy.sparse.csr_matrix(a_dense)
    
    print_sparse_blocks(extract_sparse_blocks(a))
    
def test_extract_sparse_blocks_random(n, m):
    import itertools
    
    a_dense = numpy.zeros((n, n))
    for i in xrange(n):
        if numpy.random.rand()<0.2:
            jj = min(min(i, n-1-i), 3)
            j = numpy.random.randint(-jj, jj+1)
            a_dense[i+j, i-j] = 1
    a = scipy.sparse.csr_matrix(a_dense)
    
    blocks = extract_sparse_blocks(a)
    
    print_sparse_blocks(blocks)
    display_sparse_blocks(blocks)

def test_extract_sparse_blocks_b():
    test_extract_sparse_blocks_random(100, 50)


if __name__ == '__main__':
    test_extract_sparse_blocks_a()
    test_extract_sparse_blocks_b()
