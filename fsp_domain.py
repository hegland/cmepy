import numpy

import lexarrayset

from fsp_structures import (BoundaryMatrixFactory,
                            FluxMatrixSum,
                            InteriorMatrixFactory,
                            ShiftStateTracker,
                            StateEnumeration, )

class FspDomain(object):
    def __init__(self,
                 dim,
                 reactions,
                 projections):
        
        # the dimension of the lattice that contains this domain
        self.dim = dim
        
        # enumeration of all the states in the state space
        # since the state space only ever grows, this is
        # monotonic
        self.state_enum = StateEnumeration(self.dim)
        
        self.shifts = set()
        for reaction in reactions:
            self.shifts.add(reaction.shift) 
        
        # track the sink states, ie, the reachable destination states
        # that are not included in the domain
        self.sink_states = lexarrayset.empty(self.dim)
        
        # used to track the interior and boundary source indices
        # for each shift
        self.shift_state_trackers = {}
        for shift in self.shifts:
            shift_state_tracker = ShiftStateTracker(shift, self.dim)
            self.shift_state_trackers[shift] = shift_state_tracker
        
        # for each projection, maintain an enumeration of the
        # projected sink states
        self.projections = projections
        self.proj_state_enums = {}
        for projection_name in projections:
            proj_state_enum = StateEnumeration(self.dim)
            self.proj_state_enums[projection_name] = proj_state_enum
        
        # for reach reaction, maintain an interior matrix factory
        # to produce flux matrices for the interior flux
        # (ie flow from the domain to the domain)
        self.reactions = reactions
        self.interior_matrix_factories = {}
        for reaction_name in self.reactions:
            reaction = self.reactions[reaction_name]
            shift_state_tracker = self.shift_state_trackers[reaction.shift]
            source_states = shift_state_tracker.interior
            matrix_factory = InteriorMatrixFactory(self.state_enum,
                                                   reaction.shift,
                                                   reaction.propensity)
            self.interior_matrix_factories[reaction.name] = matrix_factory        
        
    def extend(self, sigma, diff_eqs_factory):
        """
        extend the truncated domain by the states in the LexArraySet sigma.
        sigma must be disjoint with the current truncated domain.
        """
        
        # extend the state enumeration for the truncated state space
        self.state_enum.extend(sigma)
        
        sink_add = lexarrayset.empty(self.dim)
        sink_remove = lexarrayset.empty(self.dim)
        shift_updates = {}
        for shift in self.shifts:
            # extend the shift state tracker, store the resulting updates
            update = self.shift_state_trackers[shift].extend(sigma)
            shift_updates[shift] = update
            
            # accumulate the changes in boundary to the sink states
            sink_add.union_update(update.delta_plus_boundary.shift(shift))
            sink_remove.union_update(update.delta_minus_boundary.shift(shift))
        
        # apply the boundary changes to the sink states
        self.sink_states.difference_update(sink_remove)
        self.sink_states.union_update(sink_add)
        
        # completely rebuild the projection state enumerations
        # xxx todo consider optimisations maybe? (profile first)
        
        for projection_name in self.proj_state_enums:
            projection = self.projections[projection_name]
            projected_states = lexarrayset.create(projection(self.sink_states))
            proj_state_enum = self.proj_state_enums[projection_name]
            proj_state_enum.reinitialise(projected_states)
        
        # specify enumeration offsets
        # xxx todo put this as its own function / method
        # use it to generate pack / unpack functions during the
        # same pass
        
        net_offset = 0
        self.state_enum.offset = net_offset
        net_offset += self.state_enum.size
        for projection_name in self.proj_state_enums:
            proj_state_enum = self.proj_state_enums[projection_name]
            proj_state_enum.offset = net_offset
            net_offset += proj_state_enum.size
        
        matrix_size = net_offset
        
        # figure out new interior flux terms
        for reaction_name in self.reactions:
            matrix_factory = self.interior_matrix_factories[reaction_name]
            source_states = shift_updates[matrix_factory.shift].delta_interior
            interior_term = matrix_factory.produce_term(source_states,
                                                        matrix_size)
            diff_eqs_factory.update_interior_matrix(reaction_name,
                                                    interior_term)
        
        for reaction_name in self.reactions:
            pass
        
        # boundary flux :
        
        # enumerate over shifts
        # -- figure out flow out boundary term
        # -- enumerate over projections
        # -- -- figure out projection flow in boundary term
    
    def update_diff_eqs_factory(self, diff_eqs_factory):
        for reaction_name in self.reactions:
            matrix_factory = self.interior_matrix_factories[reaction_name]
        
        