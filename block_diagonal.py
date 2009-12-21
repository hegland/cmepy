import itertools
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
    def __init__(self, shape):
        assert shape is not None
        assert len(shape) == 2
        assert not (shape[0]<0 or shape[1]<0)
        self.shape = shape
        self.rows = []
        self.cols = []
        self.vals = []
    
    def add_dense_block(self,
                        block,
                        row_offset,
                        col_offset):
        
        block_rows, block_cols = numpy.indices(numpy.shape(block))
        block_rows = numpy.ravel(block_rows) + row_offset
        block_cols = numpy.ravel(block_cols) + col_offset
        block_vals = numpy.ravel(block)
        block_support = (block_vals != 0.0)
        self.rows.append(block_rows[block_support])
        self.cols.append(block_cols[block_support])
        self.vals.append(block_vals[block_support])
    
    def add_coo_block(self,
                        block,
                        row_offset,
                        col_offset):
        
        self.rows.append(block.row + row_offset)
        self.cols.append(block.col + col_offset)
        self.vals.append(block.data)
    
    def to_coo_matrix(self):
        rows = _join_arrays(self.rows)
        cols = _join_arrays(self.cols)
        vals = _join_arrays(self.vals)
        return scipy.sparse.coo_matrix((vals, (rows, cols)), self.shape)

class BlockDiagonalMatrix(object):
    def __init__(self, shape):
        assert shape is not None
        assert len(shape) == 2
        assert not (shape[0]<0 or shape[1]<0)
        self.shape = shape
        
        self.blocks = []
    
    def add_block(self, start, size, data):
        self.blocks.append((start, size, data))
    
    def add_zero_blocks(self, start, number):
        for i in xrange(number):
            data = scipy.sparse.coo_matrix((1, 1))
            self.blocks.append((start+i, 1, data))
    
    def complete_zero_blocks(self):
        block_starts = []
        block_ends = []
        for start, size, data in self.blocks:
            block_starts.append(start)
            block_ends.append(start+size-1)
        
        block_starts = numpy.array(block_starts)
        block_ends = numpy.array(block_ends)
        block_ordering = numpy.argsort(block_starts)
        block_starts = block_starts[block_ordering]
        block_ends = block_ends[block_ordering]
        
        # initial gap, if any
        if block_starts[0] > 0:
            self.add_zero_blocks(0, block_starts[0])
        # middle gaps, if any
        for prev_end, curr_start in itertools.izip(block_ends[:-1],
                                                   block_starts[1:]):
            if curr_start > prev_end+1:
                self.add_zero_blocks(prev_end+1,
                                     curr_start-prev_end-1)
        # final gap, if any
        n = min(self.shape[0], self.shape[1])
        if block_ends[-1] < n-1:
            self.add_zero_blocks(block_ends[-1]+1,
                                 n-1-block_ends[-1])
        
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
        return BlockDiagonalMatrix(a.shape)
    
    # min and max coords
    min_coord = numpy.minimum(row, col)
    max_coord = numpy.maximum(row, col)
    
    # reorder entries with respect to increasing min coord
    min_coord_order = numpy.argsort(min_coord)
    
    row = row[min_coord_order]
    col = col[min_coord_order]
    data = data[min_coord_order]
    min_coord = min_coord[min_coord_order]
    max_coord = max_coord[min_coord_order]
    
    # determine indices where blocks end
    accum_max = numpy.maximum.accumulate(max_coord)
    end_block = numpy.zeros(numpy.shape(accum_max))
    end_block[-1] = True
    end_block[:-1] = accum_max[:-1] < min_coord[1:]
    end_block_indices = numpy.nonzero(end_block)[0]
    
    # partition entries into blocks, return as BlockDiagonalMatrix
    result = BlockDiagonalMatrix(a.shape)
    block_i = 0
    for block_j in end_block_indices:
        # block_start and block_end index where this block occurs
        # along the diagonal of the original matrix a
        block_start = min_coord[block_i]
        block_end = accum_max[block_j]
        block_size = 1 + block_end - block_start
        # select the subset of the entries for this block and
        # transform into coo matrix representation
        block_slice = slice(block_i, block_j+1, None)
        coo_data = (data[block_slice],
                    (row[block_slice]-block_start,
                     col[block_slice]-block_start))
        result.add_block(block_start,
                         block_size,
                         scipy.sparse.coo_matrix(coo_data))
        block_i = block_j+1
    return result
    
    


def map(block_diagonal_matrix, f):
    result = BlockDiagonalMatrix(block_diagonal_matrix.shape)
    for block_start, block_end, block_data in block_diagonal_matrix.blocks:
        result.add_block(*f(block_start, block_end, block_data))
    return result

def expm(block_diagonal_matrix, t):
    """
    XXX TODO BUG : THIS DOESNT HANDLE ZEROS ALONG THE DIAGONAL PROPERLY
    CLEANEST WAY IS PROBABLY TO ADD THEM IN IN THE BLOCK DIAG CREATION
    ROUTINE 'from_sparse_matrix'
    """
    def _expm_block(start, end, data):
        data_dense = data.todense()*t
        exp = scipy.linalg.expm(data_dense)
        return start, end, scipy.sparse.coo_matrix(exp)
    # ensure all 'zero blocks' along the diagonal are explicitly added
    block_diagonal_matrix.complete_zero_blocks()
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
    for i, svd_block in enumerate(block_diagonal_svd.blocks):
        start, end, (u, s, vh) = svd_block
        net_s.append(s)
        block_indices = i*numpy.ones(numpy.shape(s), dtype=numpy.int)
        net_block_indices.append(block_indices)
    net_s = _join_arrays(net_s)
    net_block_indices = _join_arrays(net_block_indices)
    sort_indices = numpy.argsort(net_s)
    large_block_indices = net_block_indices[sort_indices[-k:]]
    
    # block_ks[i] = number of s values to take from the i-th svd block
    counts = numpy.bincount(large_block_indices)
    block_ks = numpy.zeros((len(block_diagonal_svd.blocks), ), dtype=numpy.int)
    block_ks[:numpy.size(counts)] = counts
    return block_ks

def to_sparse(block_diagonal_matrix):
    accumulator = Accumulator(block_diagonal_matrix.shape)
    for block_start, block_size, block_data in block_diagonal_matrix.blocks:
        accumulator.add_coo_block(block_data,
                                  block_start,
                                  block_start)
    return accumulator.to_coo_matrix()
    
def to_sparse_rank_k_approx(block_diagonal_svd, k):
    block_ks = svd_block_ks(block_diagonal_svd, k)
    block_svds = block_diagonal_svd.blocks
    f_accumulator = Accumulator(block_diagonal_svd.shape)
    e_accumulator = Accumulator(block_diagonal_svd.shape)
    for block_k, svd in itertools.izip(block_ks, block_svds):
        start, size, (u, s, vh) = svd
        if block_k == 0:
            continue
        else:           
            f_block = u[:, :block_k]
            e_block = numpy.dot(numpy.diag(s[:block_k]), vh[:block_k, :])
            f_accumulator.add_dense_block(f_block, start, start)
            e_accumulator.add_dense_block(e_block, start, start)
    f_matrix = f_accumulator.to_coo_matrix()
    e_matrix = e_accumulator.to_coo_matrix()
    return (e_matrix, f_matrix)