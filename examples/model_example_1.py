# Example 1a.
# A minimal model for reaction A -> B, where there are
# initially 10 copies of A, and 0 copies of B.
# The fields included are:
#
#    propensities    : a sequence of propensity functions, mapping reaction
#                      counts to reaction propensities, for each reaction.
#    np              : a sequence of positive integers specifying the
#                      dimensions of the reaction count state space.
#                      In this example, np contains only one element,
#                      indicating that there is only one reaction in the
#                      system, and that the corresponding reaction count
#                      takes 11 possible values: 0, 1, 2, 3, ..., 10 .
#
model_example_1a = {'propensities' : (lambda x : 1.0*(10-x), ),
                    'np' : (11, ),
                    }
# Example 1b.
# A full version of the above model, with the additional fields:
#
#    doc            : a short description of the model.
#    reactions      : sequence of reaction names, corresponding to the
#                     ordering of the propensity functions.
#    species counts : sequence of functions mapping reaction counts to
#                     species counts, for each species in the system.
#    species        : sequence of species names, corresponding to the
#                     ordering of the species count functions.
#
model_example_1b = {'doc' : 'example model for reaction A -> B',
                    'propensities' : (lambda x : 1.0*(10-x), ),
                    'reactions' : ('A->B', ),
                    'species counts' : (lambda x : 10-x,
                                        lambda x : x, ),
                    'species' : ('A', 'B', ),
                    'np' : (11, ),
                    }