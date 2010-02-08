"""
experimental cme_solver implementation
"""

import numpy

import cmepy.new_core.ode_solver as ode_solver
import cmepy.new_core.cme_matrix as cme_matrix
import cmepy.new_core.state_enum as state_enum
import cmepy.new_core.domain as domain

def create_packing_functions(domain_enum):
    """
    create_packing_functions(domain_enum) -> (pack, unpack)
    
    where
    
        pack((p, p_sink)) -> y
        unpack(y) -> (p, p_sink)
    """
    
    def pack((p, p_sink)):
        """
        pack((p, p_sink)) -> y
        
        where
        
            p : mapping from states to probability
            p_sink : float, storing probability lost from domain due to
                truncation of domain states
            y : array passed to differential equations solver
        """
        d_dense = domain_enum.pack_distribution(p)
        return numpy.concatenate((d_dense, [p_sink]))
    def unpack(y):
        """
        unpack(y) -> (p, p_sink)
        
        where
        
            p : mapping from states to probability
            p_sink : float, storing probability lost from domain due to
                truncation of domain states
            y : array passed to differential equations solver
        """
        p_sparse = domain_enum.unpack_distribution(y[:-1])
        p_sink = y[-1]
        return p_sparse, p_sink
    
    return (pack, unpack)

def create_cme_solver(model,
                      sink,
                      p_0=None,
                      time_dependencies=None,
                      domain_states=None):
    """
    create_cme_solver(model,sink[,p_0,time_dependencies,states]) -> solver
    
    returns a solver for the Chemical Master Equation of the given model.
    
    arguments:
    
        model : the CME model to solve
        
        sink : If sink is True, the solver will include a 'sink' state used
            to accumulate any probability that may flow outside the domain.
            This can be used to measure the error in the solution due to
            truncation of the domain. If sink is False, the solver will not
            include a 'sink' state, and probability will be artificially
            prevented from flowing outside of the domain.
        
        p_0 : (optional) mapping from states in the domain to probabilities,
            for the initial probability distribution. If not specified,
            and the origin of the state space can be inferred, defaults
            to all probability concentrated at the origin, otherwise, a
            ValueError will be raised.
        
        time_dependencies : (optional) mapping of time dependent coefficient
            functions keyed by subsets of reaction indices, with respect to the
            ordering of reactions determined by the order of the propensity
            functions inside the model. The propensities of the reactions
            with indices included in the subset are multiplied by the time
            dependent coefficient functions. By default, no time dependent
            coefficient functions are specified, that is, the CME has
            time-independent propensities.
        
        domain_states : (optional) array of states in the domain.
            By default, attempt to infer the domain states assuming a
            rectangular domain defined by the 'np' entry of the model, and
            optionally also the 'norigin' entry. A ValueError is raised if both
            domain_states and model['np'] are unspecified.
    """
    
    assert type(sink) is bool
    
    origin = model.get('norigin', None)
    
    # determine states in domain, then construct an enumeration of the
    # domain states
    if domain_states is None:
        if 'np' not in model:
            lament = 'if no states given, model must contain key \'np\''
            raise KeyError(lament)
        else:
            # origin is now well-defined
            if origin is None:
                origin = (0,)*len(model['np'])
            domain_states = domain.from_rect(shape = model['np'],
                                             slices = None,
                                             origin = origin)
    
    domain_enum = state_enum.create(domain_states)
    
    # determine p_0, then construct a dense representation with respect to
    # the domain enumeration
    if p_0 is None:
        if origin is None:
            lament = 'if no p_0 given, model must contain key \'norigin\''
            raise ValueError(lament)
        else:
            p_0 = {origin : 1.0}
    
    # compute reaction matrices and use them to define differential equations
    gen_matrices = cme_matrix.gen_reaction_matrices(model,
                                                    domain_enum,
                                                    sink,
                                                    cme_matrix.non_neg_states)
    reaction_matrices = list(gen_matrices)
    dy_dt = cme_matrix.create_diff_eqs(reaction_matrices,
                                       phi = time_dependencies)
    
    # construct and initialise solver
    if sink:
        sink_p_0 = 0.0
        cme_solver = ode_solver.Solver(dy_dt, y_0 = (p_0, sink_p_0))
        pack, unpack = create_packing_functions(domain_enum)
        cme_solver.set_packing(pack, unpack, transform_dy_dt = False)
    else:
        pack = domain_enum.pack_distribution
        unpack = domain_enum.unpack_distribution
        cme_solver = ode_solver.Solver(dy_dt, y_0 = p_0)
        cme_solver.set_packing(pack, unpack, transform_dy_dt = False)
    return cme_solver