import itertools
import numpy
import pylab
import enthought.mayavi.mlab as mv

def parse_solution_file(f):
    p = {}
    for line in f:
        tokens = line.split(':')
        tokens = [tok.strip() for tok in tokens]
        assert len(tokens) == 2
        mass = float(tokens[1])
        state = tokens[0]
        state = state.strip('()')
        coords = state.split(',')
        coords = tuple([int(coord.strip()) for coord in coords])
        p[coords] = mass
    return p

def plot_densities(files):
    dims = [2,3,4]
    for f in files:
        p = parse_solution_file(f)
        dense_contour_plot(p, dims)

def dense_contour_plot(p, dims):
    """
    draw 3d contour plots of 3d marginals ...
    
    NB to get mayavi threading to work, run this script via
    ipython -wthread
    """
    
    proj_mask = numpy.array([i in dims for i in xrange(5)], dtype=numpy.bool)
    cat_mask = numpy.logical_not(proj_mask)
    
    assert len(dims) == 3
    p_cat = {}
    min_coords = {}
    max_coords = {}
    for state, mass in p.iteritems():
        state = numpy.array(state, dtype=numpy.int)
        proj_state = tuple(state[proj_mask])
        cat_state = tuple(state[cat_mask])
        p_bar = p_cat.get(cat_state, {})
        p_bar[proj_state] = p_bar.get(proj_state, 0.0) + mass
        p_cat[cat_state] = p_bar
        for i, coord in enumerate(proj_state):
            if i not in min_coords:
                min_coords[i] = coord
            else:
                min_coords[i] = min(coord, min_coords[i])
            if i not in max_coords:
                max_coords[i] = coord
            else:
                max_coords[i] = max(coord, max_coords[i])
    
    proj_dims = range(len(dims))
    shape = tuple([max_coords[i] + 1 - min_coords[i] for i in proj_dims])
    
    for cat in p_cat:
        print 'category %s' % str(cat)
        
        p_bar = p_cat[cat]
        dense_p_bar = numpy.zeros(shape, dtype=numpy.float)
        print 'building dense distribution'
        for state, mass in p_bar.iteritems():
            index = tuple(numpy.array([state[i] - min_coords[i] for i in proj_dims]))
            dense_p_bar[index] = mass
        print 'plotting'
        
        mv.figure(str(cat))
        plot = mv.contour3d(dense_p_bar, contours=15, transparent=True)

def compute_marginals(p):
    marginals = {}
    for state in p:
        mass = p[state]
        for i, coord in enumerate(state):
            if i not in marginals:
                marginals[i] = {}
            marginals[i][coord] = marginals[i].get(coord, 0.0) + mass
    return marginals

def compute_expected_value(marginal):
    expected_value = 0.0
    for state in marginal:
        mass = marginal[state]
        expected_value += state*mass
    return expected_value

def compute_percentile(marginal, q):
    q = q / 100.0
    
    n = len(marginal)
    states = numpy.zeros((n, ), dtype = numpy.int)
    values = numpy.zeros((n, ), dtype = numpy.float)
    for i, (state, value) in enumerate(marginal.iteritems()):
        states[i] = state
        values[i] = value
    sort = numpy.argsort(states)
    states = states[sort]
    values = values[sort]
    
    mass = 0.0
    for i in xrange(n):
        if (mass < q) and (mass + values[i] >= q):
            if i == 0:
                state = states[i]
            else:
                gamma = (q-mass)/values[i]
                state = (1.0-gamma)*states[i-1] + gamma*states[i]
            return state
        else:
            mass += values[i]
    return states[-1]
    

def plot_marginals(files):
    def reflect(x):
        return 100-x
    percentile_pairs = [(x, reflect(x)) for x in xrange(1, 51, 1)]
    percentiles = set()
    for pair in percentile_pairs:
        percentiles.update(set(pair))
    
    percentile_values = {}
    for q in percentiles:
        percentile_values[q] = {}
    expected_values = {}
    for f in files:
        p = parse_solution_file(f)
        marginals = compute_marginals(p)
        for i in marginals:
            marginal = marginals[i]
            # expected value
            expected_value = compute_expected_value(marginal)
            ev = expected_values.get(i, [])
            ev.append(expected_value)
            expected_values[i] = ev
            # percentiles
            for q in percentiles:
                m_q = compute_percentile(marginal, q)
                qs = percentile_values[q].get(i, [])
                qs.append(m_q)
                percentile_values[q][i] = qs
    coords = marginals.keys()
    coords.sort()
    
    line_colour = 'k'
    ev_style = '-'
    q_style = '--'
    
    for i in coords:
        pylab.figure()
        pylab.plot(expected_values[i], line_colour+ev_style)
        for q_lo, q_hi in percentile_pairs:
            q_lo_data = percentile_values[q_lo][i]
            q_hi_data = percentile_values[q_hi][i]
            n = len(q_lo_data)
            x_lo = range(n)
            x_hi = range(n)
            x_lo.reverse()
            x = x_lo + x_hi
            y_lo = list(q_lo_data)
            y_hi = list(q_hi_data)
            y_lo.reverse()
            y = y_lo + y_hi
            pylab.fill(x, y, line_colour+q_style, facecolor='blue', alpha=0.05)
        pylab.savefig('marginal_%04d.png' % i)
        pylab.close()

def main():
    def solution_file(i):
        file_name = 'solution_%04d.txt' % i
        return file(file_name, 'r')
    files = [solution_file(i) for i in xrange(10)]
    #plot_marginals(files)
    plot_densities(files)
    for f in files:
        f.close()

if __name__ == '__main__':
    main()
