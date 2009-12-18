import numpy
import scipy.sparse

def _join_arrays(arrays):
    net_size = 0
    for a in arrays:
        net_size += numpy.size(a)
    
    result = numpy.zeros((net_size, ), dtype = arrays[0].dtype)
    i = 0
    for a in arrays:
        size = numpy.size(a)
        result[i:i+size] = a
        i += size
    return result

class Accumulator(object):
    def __init__(self):
        self.rows = []
        self.cols = []
        self.vals = []
        self.size = 0
    
    def add_dense_block(self,
                        block,
                        row_offset,
                        col_offset,
                        block_size):
       
        assert (block_size >= 0)
        if block_size == 0:
            return
        
        block_rows, block_cols = numpy.indices(numpy.shape(block))
        block_rows = numpy.ravel(block_rows) + row_offset
        block_cols = numpy.ravel(block_cols) + col_offset
        block_vals = numpy.ravel(block)
        block_support = numpy.ravel(block != 0.0)
        self.rows.append(block_rows[block_support])
        self.cols.append(block_cols[block_support])
        self.vals.append(block_vals[block_support])
        self.size += block_size
    
    def add_coo_block(self,
                        block,
                        row_offset,
                        col_offset,
                        block_size):
        
        assert (block_size >= 0)
        if block_size == 0:
            return
        
        self.rows.append(block.row + row_offset)
        self.cols.append(block.col + col_offset)
        self.vals.append(block.data)
        self.size += block_size
    
    def to_coo_matrix(self):
        rows = _join_arrays(self.rows)
        cols = _join_arrays(self.cols)
        vals = _join_arrays(self.vals)
        shape = (self.size, self.size)
        return scipy.sparse.coo_matrix((vals, (rows, cols)), shape)

def from_sparse_matrix(a):
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


def map(block_diagonal_matrix, f):
    result = []
    for block_start, block_end, block_data in block_diagonal_matrix:
        result.append(f(block_start, block_end, block_data))
    return result

def expm(block_diagonal_matrix, t):
    """
    XXX TODO BUG : THIS DOESNT HANDLE ZEROS ALONG THE DIAGONAL PROPERLY
    CLEANEST WAY IS PROBABLY TO ADD THEM IN IN THE BLOCK DIAG CREATION
    ROUTINE 'from_sparse_matrix'
    """
    def _expm_block(start, end, data):
        data_dense = data.todense()*t
        exp = scipy.linalg.expm(data_dense, q = 7)
        return start, end, scipy.sparse.coo_matrix(exp)
    return map(block_diagonal_matrix, _expm_block)

def block_svd(block_diagonal_matrix):
    def _svd_block(start, end, data):
        data_dense = data.todense()
        u, s, vh = scipy.linalg.svd(data_dense)
        return start, end, (u, s, vh)
    return map(block_diagonal_matrix, _svd_block)

def svd_block_ks(block_diagonal_svd, k):
    net_s = []
    net_block_indices = []
    
    # figure out how many s values to take from each block,
    # given that we are only interested in the k largest s values
    # over all blocks.
    for i, svd_block in enumerate(block_diagonal_svd):
        start, end, (u, s, vh) = svd_block
        net_s.append(s)
        block_indices = i*numpy.ones(numpy.shape(s), dtype=numpy.int)
        net_block_indices.append(block_indices)
    net_s = _join_arrays(net_s)
    net_block_indices = join_arrays(net_block_indices)
    sort_indices = numpy.argsort(net_s)
    large_block_indices = net_block_indices[sort_indices[-k:]]
    
    # block_ks[i] = number of s values to take from the i-th svd block
    counts = numpy.bincount(large_block_indices)
    block_ks = numpy.zeros((len(block_diagonal_svd), ), dtype=numpy.int)
    block_ks[:numpy.size(counts)] = counts
    return block_ks

def to_sparse(block_diagonal_matrix):
    accumulator = Accumulator()
    for block_start, block_size, block_data in block_diagonal_matrix:
        accumulator.add_coo_block(block_data,
                                  block_start,
                                  block_start,
                                  block_size)
    return accumulator.to_coo_matrix()
    
def to_sparse_rank_k_approx(block_diagonal_svd, k):
    block_ks = svd_block_ks(block_diagonal_svd, k)
    f_accumulator = Accumulator()
    e_accumulator = Accumulator()
    for block_k, svd_block in itertools.izip(block_ks, block_diagonal_svd):
        start, size, (u, s, vh) = svd_block
        if block_k == 0:
            continue
        else:           
            f_block = u[:, :block_k]
            e_block = numpy.dot(numpy.diag(s[:block_k]), vh[:block_k, :])
            f_accumulator.add_dense_block(f_block, start, start, size)
            e_accumulator.add_dense_block(e_block, start, start, size)
    f_matrix = f_accumulator.to_coo_matrix()
    e_matrix = e_accumulator.to_coo_matrix()
    return (e_matrix, f_matrix)