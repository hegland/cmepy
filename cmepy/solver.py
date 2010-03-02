"""
creates solvers for the chemical master equation
"""

import numpy
from cmepy import cme_matrix, domain, ode_solver, state_enum
from cmepy import model as mdl

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

def create(model,
           sink,
           p_0=None,
           t_0=None,
           time_dependencies=None,
           domain_states=None):
    """
    create(model,sink[,p_0,time_dependencies,states]) -> solver
    
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
            and the initial state of the state space is given by the model,
            defaults to all probability concentrated at the initial state,
            otherwise, a ValueError will be raised.
        
        t_0 : (optional) initial time, defaults to 0.0
        
        time_dependencies : (optional) By default, reaction propensities are
            time independent. If specified, time_dependencies must be of the
            form { s_1 : phi_1, ..., s_n : phi_n }, where each (s_j, phi_j)
            item satisifes :
            
                s_j : set of reaction indices
                phi_j : phi_j(t) -> time dependent coefficient 
            
            The propensities of the reactions with indicies contained in s_j
            will all be multiplied by the coefficient phi_j(t), at time t.
            Reactions are indexed according to the ordering of the propensities
            in the model.
            
            The reaction index sets s_j must be *disjoint*. It is not necessary
            for the union of the s_j to include all the reaction indices.
            If a reaction's index is not contained in any s_j then the reaction
            is treated as time-independent. 
        
        mapping of time dependent coefficient
            functions keyed by subsets of reaction indices, with respect to the
            ordering of reactions determined by the order of the propensity
            functions inside the model. The propensities of the reactions
            with indices included in each subset are multiplied by the time
            dependent coefficient function. By default, no time dependent
            coefficient functions are specified, that is, the CME has
            time-independent propensities.
        
        domain_states : (optional) array of states in the domain.
            By default, generate the rectangular lattice of states defined by
            the 'shape' entry of the model. A ValueError is raised if both
            domain_states and 'shape' are unspecified.
    """
    
    mdl.validate_model(model)
    
    
    # determine states in domain, then construct an enumeration of the
    # domain states
    if domain_states is None:
        if mdl.SHAPE not in model:
            lament = 'if no states given, model must contain key \'%s\''
            raise KeyError(lament % mdl.SHAPE)
        else:
            domain_states = domain.from_rect(shape = model.shape)
    
    domain_enum = state_enum.create(domain_states)
    
    # determine p_0, then construct a dense representation with respect to
    # the domain enumeration
    initial_state = model.get(mdl.INITIAL_STATE, None)
    if p_0 is None:
        if initial_state is None:
            lament = 'if no p_0 given, model must contain key \'%s\''
            raise ValueError(lament % mdl.INITIAL_STATE)
        else:
            p_0 = {initial_state : 1.0}
    
    if t_0 is None:
        t_0 = 0.0
    
    member_flags = domain_enum.contains(domain.from_iter(p_0))
    if not numpy.logical_and.reduce(member_flags):
        raise ValueError('support of p_0 is not a subset of domain_states')
    
    # compute reaction matrices and use them to define differential equations
    gen_matrices = cme_matrix.gen_reaction_matrices(
        model,
        domain_enum,
        sink,
        cme_matrix.non_neg_states
    )
    reaction_matrices = list(gen_matrices)
    dy_dt = cme_matrix.create_diff_eqs(
        reaction_matrices,
        phi = time_dependencies
    )
    
    # construct and initialise solver
    if sink:
        sink_p_0 = 0.0
        cme_solver = ode_solver.Solver(
            dy_dt,
            y_0 = (p_0, sink_p_0),
            t_0 = t_0
        )
        pack, unpack = create_packing_functions(domain_enum)
        cme_solver.set_packing(
            pack,
            unpack,
            transform_dy_dt = False
        )
    else:
        pack = domain_enum.pack_distribution
        unpack = domain_enum.unpack_distribution
        cme_solver = ode_solver.Solver(
            dy_dt,
            y_0 = p_0,
            t_0 = t_0
        )
        cme_solver.set_packing(
            pack,
            unpack,
            transform_dy_dt = False
        )
    return cme_solver
