import itertools
import numpy
import fsp

import cmepy.new_core.ode_solver as ode_solver

def vectorise_states(states):
    vectstates = [numpy.array(coord) for coord in itertools.izip(*list(states))]
    vectstates = [numpy.ravel(i) for i in vectstates]
    return vectstates

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



class FspSolver(object):
    def __init__(self,
                 model,
                 p_0,
                 time_dependencies = None):
        
        self.model = model
        self.p_0 = p_0
        self.time_dependencies = time_dependencies
        
        self.error_projections = None
        self.error_projection_dims = None
        self.delta_t = None
        self.num_steps = None
        self.num_mini_steps = None
        self.tolerance = None
        self.gamma = None
        self.rho = None

    def set_error_projections(self, error_projections, error_projection_dims):
        
        self.error_projections = error_projections
        self.error_projection_dims = error_projection_dims
    
    def set_fsp_parameters(self,
                           delta_t,
                           num_steps,
                           num_mini_steps,
                           tolerance,
                           gamma = 1.5,
                           rho = 0.9):
        
        self.delta_t = float(delta_t)
        assert (self.delta_t > 0.0), 'time step size delta_t must be positive'
        self.num_steps = int(num_steps)
        assert (self.num_steps >= 0), 'num_steps must be non-negative'
        self.num_mini_steps = int(num_mini_steps)
        assert (self.num_mini_steps > 0), 'num_mini_steps must be positive'
        self.tolerance = float(tolerance)
        assert (0.0 < self.tolerance < 1.0), 'tolerance out of range (0, 1)'
        self.gamma = float(gamma)
        assert (1.0 < self.gamma), 'gamma out of range (1, +infty)'
        self.rho = float(rho)
        assert (0.0 < self.rho < 1.0), 'rho out of range (0, 1)'

    def approx_step_truncated_cme(self, p_0, t_0, epsilon):
        print '[approx_step_truncated_cme]'
        print '\tt_0 : %g' % t_0
        print '\tepsilon : %g' % epsilon
        states = support(p_0)
        success = False
        while not success:
            success, p, expand_states = self.step_cme(p_0,
                                                      t_0,
                                                      epsilon,
                                                      states)
            states.update(expand_states)
        return p
    
    def solve(self, solution_file_prefix):
        t, p = 0.0, self.p_0
        epsilon = self.tolerance * (self.gamma**(-self.num_steps))
        
        def write_solution(p, step):
            file_name = '%s_solution_%04d.txt' % (solution_file_prefix, step)
            f = file(file_name, 'w')
            write_distribution(p, f)
            f.close()
        
        write_solution(p, 0)
        
        for step in xrange(self.num_steps):
            print '[solve_fsp]'
            print '\tstep %d, gamma %g, epsilon %g' % (step, self.gamma, epsilon)
            
            # step the solution forwards using the fsp within tolerance
            epsilon_step = self.rho*(self.gamma - 1.0)*epsilon
            p = self.approx_step_truncated_cme(p, t, epsilon_step)
            
            # compress the solution with tolerance
            epsilon_compress = (1.0 - self.rho)*(self.gamma - 1.0)*epsilon
            p = compress(p, epsilon_compress)
            
            epsilon = self.gamma*epsilon
            t += self.delta_t
            
            write_solution(p, step+1)
    
    def create_boundary_states(self, state_index_map):
        boundary_states = set()
        
        for offset_vector in self.model['offset_vectors']:
            offset_vector = numpy.array(offset_vector)[:, numpy.newaxis]
            dest_states = state_index_map.vectstates + offset_vector
            
            i_mask = numpy.logical_not(state_index_map.vectstates_in_map(dest_states))
            j_mask = numpy.logical_and.reduce(dest_states[:, i_mask]>=0, axis=0)
            
            legal_boundary_states = dest_states[:, i_mask][:, j_mask]
            
            for state in itertools.izip(*legal_boundary_states):
                boundary_states.add(state)
            
        return boundary_states
            
    def create_error_states(self,
                            state_index_map_offset,
                            state_index_map):
        
        boundary_states = self.create_boundary_states(state_index_map)
        boundary_vectstates = vectorise_states(boundary_states)
        boundary_states_list = list(boundary_states)
        
        error_state_index_maps = []
        error_states_to_dest_states = []
        
        for error_projection in self.error_projections:
            #print '\tcreating error states for error projection'
            error_states = {}
            
            proj_states = error_projection(*boundary_vectstates)
            
            proj_states_list = zip(*proj_states)
            
            assert len(proj_states_list) == len(boundary_states_list)
            
            for proj_state, dest_state in itertools.izip(proj_states_list,
                                                         boundary_states_list):
                if proj_state not in error_states:
                    error_states[proj_state] = set()
                error_states[proj_state].add(dest_state)
            
            error_states_to_dest_states.append(error_states)
            num_error_states = len(error_states)
            
            error_vectstates = vectorise_states(error_states.keys())
            
            error_state_index_map = fsp.StateIndexMap(num_error_states,
                                                      error_vectstates,
                                                      state_index_map_offset)
            error_state_index_maps.append(error_state_index_map)
            state_index_map_offset += num_error_states
        
        return (error_state_index_maps,
                error_states_to_dest_states,
                state_index_map_offset)
    
    def step_cme(self, p_0_sparse, t_0, error_tol, states):
            
        #print 'creating vectstates'
        vectstates = vectorise_states(states)
        
        #print 'creating state index map'
        state_index_map_offset = 0
        state_index_map = fsp.StateIndexMap(len(vectstates),
                                            vectstates,
                                            state_index_map_offset)
        state_index_map_offset += state_index_map.size
        
        #print 'creating error states'
        
        temp = self.create_error_states(state_index_map_offset,
                                        state_index_map)
        
        (error_state_index_maps,
         error_states_to_dest_states,
         state_index_map_offset) = temp
            
        net_state_index_map_size = state_index_map.size
        for error_state_index_map in error_state_index_maps:
            net_state_index_map_size += error_state_index_map.size
        
        error_trackers = []
        for proj, map in itertools.izip(self.error_projections,
                                        error_state_index_maps):
            error_trackers.append((proj, map))
    
        flux_matrices = fsp.create_flux_matrices(self.model,
                                                 state_index_map,
                                                 error_trackers,
                                                 self.error_projection_dims)
    
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
        
        dy_dt = fsp.create_diff_eqs(net_state_index_map_size,
                                    flux_matrices,
                                    self.time_dependencies)
        
        
        p_0 = state_index_map.pack_distribution(p_0_sparse)
        initial_mass = numpy.add.reduce(p_0)

        solver = ode_solver.Solver(dy_dt, (p_0, None), t_0)
        solver.set_packing(pack, unpack, transform_dy_dt=False)
        
        
        # from memory this is the accuracy of the solver
        solver_error = self.num_mini_steps*1.0e-7
        assert error_tol > solver_error, 'error tol %g is larger than solver error %g' % (error_tol, solver_error)
        
        t_final = t_0 + self.delta_t
        time_steps = numpy.linspace(t_0, t_final, self.num_mini_steps+1)
        
        print '\t** omega states %d ; all states %d' % (len(states), net_state_index_map_size)
        print '\t   error_tol %g' % error_tol
        print '\t   initial_mass %g' % initial_mass
        
        early_cutoff_factor = 2.0
        
        for t in time_steps[1:]:
            solver.step(t)
            p, p_errors = solver.y
            
            interior_mass = numpy.add.reduce(numpy.ravel(p))
            trunc_error = initial_mass - interior_mass
            error = trunc_error + solver_error
            
            if error > error_tol*early_cutoff_factor:
                break
            
        if error > error_tol:
            print '\t-- t = %g, error exceeds tolerance' % t
            print '\t   trunc error - tolerance : %g' % (error - error_tol)
            # measure size of support of p
            ps = numpy.argsort(p)
            psa = numpy.add.accumulate(p[ps])
            sig = psa > error_tol
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
            p_sparse = state_index_map.unpack_distribution(p)
            return True, p_sparse, set()

def test_gene_toggle():
    import munsky_khammash_gene_toggle_08
    model = munsky_khammash_gene_toggle_08.create_model(0, 0)
    p_0 = {(0, 0) : 1.0}
    fsp_solver = FspSolver(model, p_0)
    # specify state projections to use for error tracking
    error_projections = (lambda s_1, s_2 : (s_1, s_2), )
    error_projection_dims = 2
    fsp_solver.set_error_projections(error_projections, error_projection_dims)
    # specify parameters
    fsp_solver.set_fsp_parameters(delta_t = 0.2,
                                  num_steps = 10,
                                  num_mini_steps = 10,
                                  tolerance = 1.0e-2,
                                  gamma = 1.5,
                                  rho = 0.9)
    # solve the model, writing solutions to file
    fsp_solver.solve('gene_toggle')

def test_transcription_regulation():
    import goutsias_transcription_regulation_2 as tr
    dna_count = 2
    model = tr.create_model(dna_count)
    p_0 = {(2, 0, 0, 2, 6) : 1.0}
    time_dependencies = tr.create_time_dependencies()
    fsp_solver = FspSolver(model, p_0, time_dependencies)
    # specify state projections to use for error tracking
    error_projection_dims = 4
    def make_proj((i, j)):
        return lambda *x : (x[0], x[1], x[i], x[j])
    proj_args = [(2, 3), (2, 4), (3, 4)]
    error_projections = [make_proj(args) for args in proj_args]
    fsp_solver.set_error_projections(error_projections, error_projection_dims)
    # specify parameters
    fsp_solver.set_fsp_parameters(delta_t = 10.0,
                                  num_steps = 10,
                                  num_mini_steps = 20,
                                  tolerance = 1.0e-2,
                                  gamma = 1.5,
                                  rho = 0.9)
    # solve the model, writing solutions to file
    fsp_solver.solve('transcription_regulation')

def profile_main():
    import cProfile, pstats
    profile_file = 'fsp_tests.profile'
    cProfile.run('test_transcription_regulation()', profile_file)
    stats = pstats.Stats(profile_file)
    stats.sort_stats('cumulative').print_stats(30)

if __name__ == '__main__':
    profile_main()