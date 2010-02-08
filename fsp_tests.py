import itertools

import numpy

import fsp
import goutsias_transcription_regulation as tr
import cmepy.new_core.ode_solver as ode_solver
import cmepy.new_core.cme_solver as cme_solver
import cmepy.new_core.recorder as cme_recorder


# xxx todo add serialisation of solutions to file
# easiest way is probably to timestamp the file name
# then dump a list of state : nonzero prob pairs

# xxx todo think about expansion states
# is there some monotonic thing going on
# eg if x and y are error states containing mass a and b
# then if x is added to the state space the error accumulated inside y
# should still be at least b ... right ?!
# this would allow an obvious optimisation to the expansion approach
# as long as the error per state estimate is decent enough

# idea for error tracking / expansion heuristic
# 
# VANILLA FSP EXPAND APPROACH : use 1 error state for all error
#    - gives no information *where* error is occurring
#    - cheap
# MY EXPAND APPROACH : use unique error state for each possible error
#    - gives maximal information where error is occurring
#    - not at all cheap
# a different expand approach :
#    - define a number of categories based on the value of a
#      single state coord and a reaction offset
#      --- eg (COORD A, VALUE B, OFFSET C)
#    accumulate error in error states for these categories
#    .... this is equivalent to using original approach
#         over all 1-variable marginal distributions
#         from the full distribution
#    .... obvious generalisation : use m-variable marginal
#         distributions for some 1 <= m <= dims
#         where dims is the dimension of the state space ...
#
#

def write_distribution(p, f):
    for state, mass in p.iteritems():
        line = '%s : %g\n' % (str(state), mass)
        f.write(line)

def support(p):
    states = set()
    for state in p:
        states.add(state)
    return states

def compress(p, epsilon):
    print '[compressing]'
    print '\tepsilon : %g' % epsilon
    print '\tpre-compression size : %d' % len(p)
    states = p.keys()
    n = len(states)
    mass = numpy.array(p.values())
    sorted = numpy.argsort(mass)
    acc_sorted_mass = numpy.add.accumulate(mass[sorted])
    p_compressed = {}
    for i in xrange(n):
        # ignore the smallest states with
        # probability summing to less than epsilon
        if acc_sorted_mass[i] < epsilon:
            continue
        state = states[sorted[i]]
        p_compressed[state] = p[state]
    print '\tpost-compression size : %d' % len(p_compressed)
    return p_compressed

def approx_step_truncated_cme(p_0, t_0, delta_t, epsilon, num_mini_steps):
    print '[approx_step_truncated_cme]'
    print '\tt_0 : %g' % t_0
    print '\tepsilon : %g' % epsilon
    states = support(p_0)
    success = False
    while not success:
        success, p, expand_states = step_cme(p_0, t_0, delta_t, epsilon, states, num_mini_steps)
        states.update(expand_states)
    return p
    
def solve_fsp(p_0, tolerance, gamma, delta_t, num_steps, rho, num_mini_steps):
    t, p = 0.0, p_0
    epsilon = tolerance * (gamma**(-num_steps))
    for step in xrange(num_steps):
        print '[solve_fsp]'
        print '\tstep %d, gamma %g, epsilon %g' % (step, gamma, epsilon)
        epsilon_step = rho*(gamma - 1.0)*epsilon
        epsilon_compress = (1.0 - rho)*(gamma - 1.0)*epsilon
        p = approx_step_truncated_cme(p, t, delta_t, epsilon_step, num_mini_steps)
        t += delta_t
        p = compress(p, epsilon_compress)
        epsilon = gamma*epsilon
        
        file_name = 'solution_%04d.txt' % step
        f = file(file_name, 'w')
        write_distribution(p, f)
        f.close()
    return p

def main():   
    p_0 = {(2, 0, 0, 2, 6) : 1.0}
    tolerance = 1.0e-2
    gamma = 1.5 # let e_n = gamma*e_(n-1)
    delta_t = 10.0
    num_steps = 10
    num_mini_steps = 21
    rho = 0.9
    p = solve_fsp(p_0, tolerance, gamma, delta_t, num_steps, rho, num_mini_steps)
    
def step_cme(p_0_sparse, t_0, delta_t, error_tol, states, num_mini_steps):
    
    dna_count = 2
    rna_max = 15
    m_max = 15
    d_max = 15
    model = tr.create_transcription_regulation_model(dna_count)
    np = (dna_count+1, )*2 + (rna_max+1, m_max+1, d_max+1)
    model['np']= np
        
    #print 'creating vectstates'
    vectstates = [numpy.array(coord) for coord in itertools.izip(*list(states))]
    vectstates = [numpy.ravel(i) for i in vectstates]
    
    #print 'creating state index map'
    state_index_map_offset = 0
    state_index_map = fsp.StateIndexMap(len(vectstates),
                                        vectstates,
                                        state_index_map_offset)
    state_index_map_offset += state_index_map.size
    
    
    projection_dim = 4
    def make_proj((i, j)):
        #return lambda *x : (x[0], x[1], x[i], x[j])
        return lambda *x : (x[0], x[1], x[i], x[j])
    proj_args = [(2, 3), (2, 4), (3, 4)]
    error_projections = [make_proj(args) for args in proj_args]
    
    #print 'creating error states'
    
    error_state_index_maps = []
    error_states_to_dest_states = []
    for error_projection in error_projections:
        #print '\tcreating error states for error projection'
        error_states = {}
        for offset_vector in model['offset_vectors']:
            offset_vector = numpy.asarray(offset_vector)[:, numpy.newaxis]
            dest_states = vectstates + offset_vector
            for dest_state in itertools.izip(*dest_states):
                if dest_state not in states:
                    # enforce non negative species counts
                    if any(itertools.imap(lambda x : x<0, dest_state)):
                        continue
                    proj_state = error_projection(*dest_state)
                    if proj_state not in error_states:
                        error_states[proj_state] = set()
                    error_states[proj_state].add(dest_state)
        
        error_states_to_dest_states.append(error_states)
        num_error_states = len(error_states)
        
        #print '\tfound %d projected error states' % num_error_states
        #print '\tvectorising projected error states'
        error_vectstates = [[] for i in xrange(projection_dim)]
        for state in error_states:
            for i, coord in enumerate(state):
                error_vectstates[i].append(coord)
        #print '\tcreating error state index map'
        error_state_index_map = fsp.StateIndexMap(num_error_states,
                                                  error_vectstates,
                                                  state_index_map_offset)
        error_state_index_maps.append(error_state_index_map)
        state_index_map_offset += num_error_states
        
    net_state_index_map_size = state_index_map.size
    for error_state_index_map in error_state_index_maps:
        net_state_index_map_size += error_state_index_map.size
    
    error_trackers = []
    for proj, map in itertools.izip(error_projections,
                                    error_state_index_maps):
        error_trackers.append((proj, map))

    flux_matrices = fsp.create_flux_matrices(model,
                                             state_index_map,
                                             error_trackers)

    # define pack / unpack as maps between dense np shaped vector p
    def pack((p, p_errors)):
        y = numpy.zeros(net_state_index_map_size)
        y[state_index_map.indices] = p
        if p_errors is not None:
            offset = state_index_map.size
            for p_error, error_state_index_map in itertools.izip(p_errors, error_state_index_maps):
                y[offset:offset+error_state_index_map.size] = p_error[:]
                offset += error_state_index_map.size
        return y
    
    def unpack(y):
        p = numpy.zeros((state_index_map.size, ))
        p[:] = y[state_index_map.indices]
        p_errors = []
        offset = state_index_map.size
        for error_state_index_map in error_state_index_maps:
            p_error = numpy.zeros((error_state_index_map.size, ))
            p_error[:] = y[offset:offset+error_state_index_map.size]
            p_errors.append(p_error)
            offset += error_state_index_map.size
        return p, p_errors
    
    time_dependencies = tr.create_transcription_regulation_time_dependencies()
    dy_dt = fsp.create_diff_eqs(net_state_index_map_size,
                                flux_matrices,
                                time_dependencies)
    
    initial_mass = 0.0
    p_0 = numpy.zeros((state_index_map.size, ))
    for state in p_0_sparse:
        index = state_index_map.state_to_index[state]
        mass = p_0_sparse[state]
        p_0[index] = mass
        initial_mass += mass
    
    solver = ode_solver.Solver(dy_dt, (p_0, None), t_0)
    solver.set_packing(pack,
                       unpack,
                       transform_dy_dt=False)
    
    
    # from memory this is the accuracy of the solver
    solver_error = 1.0e-7
    assert error_tol > solver_error, 'error tol %g is larger than solver error %g' % (error_tol, solver_error)
    # left over error tolerance for state space truncation
    
    t_final = t_0 + delta_t
    time_steps = numpy.linspace(t_0, t_final, num_mini_steps)
    
    
    print '\t** omega states %d ; all states %d' % (len(states), net_state_index_map_size)
    print '\t   error_tol %g' % error_tol
    print '\t   initial_mass %g' % initial_mass
    
    for t in time_steps[1:]:
        solver.step(t)
        p, p_errors = solver.y
        
        interior_mass = numpy.add.reduce(numpy.ravel(p))
        trunc_error = initial_mass - interior_mass
        error = trunc_error + solver_error
    if error > error_tol:
        print '\t-- t = %g, error exceeds tolerance' % t
        print '\t   trunc error - tolerance : %g' % (error - error_tol)
        # measure size of support of p
        ps = numpy.argsort(p)
        psa = numpy.add.accumulate(p[ps])
        sig = psa>error_tol
        num_sig = numpy.add.reduce(sig)
        print '\t   efficiency estimate : %g' % (float(num_sig)/len(states))
        print ''
        # figure out what the extended state space should be
        candidates = {}
        for i, p_error in enumerate(p_errors):
            for j, err in enumerate(p_error):
                proj_state = error_state_index_maps[i].index_to_state[j + error_state_index_maps[i].origin]
                pre_image_size =  len(error_states_to_dest_states[i][proj_state])
                # weight err by pre image size
                weighted_err = err / pre_image_size
                key = (weighted_err, err)
                if key not in candidates:
                    candidates[key] = set()
                candidates[key].add((i, proj_state))
        errors = candidates.keys()
        errors.sort()
        errors.reverse()
        
        extension_states = set(states)
        
        stop = False
        
        extended_error = 0.0
        extended_error_limit = 0.9*error
        
        for err_key in errors:
            if stop:
                break
            weighted_err, err = err_key
            for proj_index, proj_state in candidates[err_key]:
                proj_pre_image = error_states_to_dest_states[proj_index][proj_state]
                extension_states.update(proj_pre_image)
                extended_error += err
                if extended_error >= extended_error_limit:
                    stop = True
                    break
        extension_states = extension_states.difference(states)
        return False, None, extension_states
    else:
        print ''
        p_dict = {}
        for i in state_index_map.index_to_state:
            state = state_index_map.index_to_state[i]
            p_dict[state] = p[i]
        return True, p_dict, set()

if __name__ == '__main__':
    import cProfile, pstats
    profile_file = 'fsp_tests.profile'
    cProfile.run('main()', profile_file)
    stats = pstats.Stats(profile_file)
    stats.sort_stats('cumulative').print_stats(30)