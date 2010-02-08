import numpy

import cmepy.new_core.ode_solver as ode_solver
from cmepy.core.core_cme import process_offset_vectors, gen_reaction_actions


def create_flux_data(model):    
    np = model['np']
    num_states = numpy.prod(np)
    propensities = model['propensities']
    norigin = model.get('norigin', None)
    offset_vectors = model.get('offset_vectors', None)
    offset_vectors = process_offset_vectors(np,
                                            propensities,
                                            offset_vectors) 
    reaction_actions = gen_reaction_actions(np,
                                            propensities,
                                            offset_vectors,
                                            norigin)
    return reaction_actions

def default_flux_evaluator(reaction_index,
                           t,
                           p,
                           coefficients,
                           source_indices,
                           dest_indices):
    return coefficients*p[source_indices]

def create_time_dependent_flux_evaluator(time_dependencies):
    def time_dependent_flux_evaluator(reaction_index,
                           t,
                           p,
                           coefficients,
                           source_indices,
                           dest_indices):
        phi = time_dependencies[reaction_index](t)
        return coefficients*p[source_indices]*phi
    return time_dependent_flux_evaluator

def create_diff_eqs(flux_data, flux_evaluator=None):
    
    if flux_evaluator is None:
        flux_evaluator = default_flux_evaluator
    
    def diff_eqs(t, p):
        """
        this routine defines our mapping diff_eqs(t,p) ---> dp/dt(t,p)
        which is later called by the ode solver
        """
        
        inflated_shape = numpy.shape(p)
        
        # the net flux will be accumulated inside the array p_dot
        p_dot = numpy.zeros(inflated_shape)
        
        for reaction_index, data in enumerate(flux_data):
            
            coefficients, source_indices, dest_indices = data
            
            # compute flux for this reaction
            flux = flux_evaluator(reaction_index,
                                  t,
                                  p,
                                  coefficients,
                                  source_indices,
                                  dest_indices)
            
            # flux flows from source to destination
            p_dot[source_indices] -= flux
            p_dot[dest_indices] += flux
        
        # return the net flux
        return p_dot
    
    return diff_eqs

def create_default_p_0(model):
    inflated_shape = model['np']
    flattened_shape = (numpy.product(inflated_shape), )
    default_p_0 = numpy.zeros(flattened_shape)
    default_p_0[0] = 1.0
    return numpy.reshape(default_p_0, inflated_shape)

def create_packing_functions(model):
    inflated_shape = model['np']
    flattened_shape = (numpy.product(inflated_shape), )
    
    def pack(p):
        return numpy.reshape(p, flattened_shape)
    
    def unpack(y):
        return numpy.reshape(y, inflated_shape)
    
    return (pack, unpack)

def create_cme_solver(model):
    # construct dependencies
    pack, unpack = create_packing_functions(model)
    default_p_0 = create_default_p_0(model)
    flux_data = create_flux_data(model)
    dp_dt = create_diff_eqs(flux_data)
    # then initialise and return solver
    cme_solver = ode_solver.Solver(dp_dt, default_p_0)
    cme_solver.set_packing(pack, unpack)
    return cme_solver