# Example 2a.
# A minimal model for system of reactions A -> B, B -> C
# initially 10 copies of A, and 0 copies of both B & C.
#
model_example_2a = {'propensities' : (lambda x1, x2 : 1.0*(10-x1),
                                      lambda x1, x2 : 1.0*(x1-x2),),
                    'np' : (11, 11),
                    }
# Example 2b.
# A full version of the model 2a.
#
model_example_2b = {'doc' : 'example model for reactions A -> B, B -> C',
                    'propensities' : (lambda x1, x2 : 1.0*(10-x1),
                                      lambda x1, x2 : 1.0*(x1-x2), ),
                    'reactions' : ('A->B', 'B->C', ),
                    'species counts' : (lambda x1, x2 : 10-x1,
                                        lambda x1, x2 : x1-x2,
                                        lambda x1, x2 : x2, ),
                    'species' : ('A', 'B', 'C', ),
                    'np' : (11, 11, ),
                    }

# Example 2c.
# An alternate version of model 2b, this time using python 'star magic'.
#
model_example_2c = {'doc' : 'example model for reactions A -> B, B -> C',
                    'propensities' : (lambda *x : 1.0*(10-x[0]),
                                      lambda *x : 1.0*(x[0]-x[1]), ),
                    'reactions' : ('A->B', 'B->C'),
                    'species counts' : (lambda *x : 10-x[0],
                                        lambda *x : x[0]-x[1],
                                        lambda *x : x[1], ),
                    'species' : ('A', 'B', 'C', ),
                    'np' : (11, 11, ),
                    }

# Example 2d.
# Yet another variation of 2b, defining the propensities in terms of
# the species counts.
#

def generate_model_example_2d(initial_count_a,
                              initial_count_b,
                              initial_count_c):
    """
    Generates and returns a complete model for the system of reactions
    A->B, B->C, using the specified initial counts of A, B & C.
    """
    
    # define maps from reaction counts to species counts first
    species_counts = (lambda *x : initial_count_a - x[0],
                      lambda *x : initial_count_b + x[0] - x[1],
                      lambda *x : initial_count_c + x[1], )
    
    # then define the reaction propensities in terms of the species counts
    propensities = (lambda *x : 1.0*species_counts[0](*x),
                    lambda *x : 1.0*species_counts[1](*x), )
    
    model = {'doc' : 'example model for reactions A -> B, B -> C',
             'propensities' : propensities,
             'reactions' : ('A->B', 'B->C', ),
             'species counts' : species_counts,
             'species' : ('A', 'B', 'C', ),
             'np' : (initial_count_a + 1,
                     initial_count_a + initial_count_b + 1, ),
             }
    
    return model

model_example_2d = generate_model_example_2d(10, 0, 0)


def test_models():
    """
    Runs a unit test to verify that the models defined in the module all
    lead to the same solution, at t = 10.0 seconds, when solved via the
    CmeSolver.
    """
    
    import unittest
    from test import test_support

    class ModelTest(unittest.TestCase):
        def setUp(self):
            self.models = [model_example_2a,
                           model_example_2b,
                           model_example_2c,
                           model_example_2d,
                           ]
        
        def testModelEquivalence(self):
            import numpy
            import cmepy.solver
            from numpy.testing.utils import assert_array_almost_equal
            
            time_steps = numpy.linspace(0.0, 10.0, 11)
            prev_soln = None
            for model in self.models:
                solver = cmepy.solver.CmeSolver(model)
                for t in time_steps:
                    solver.step(t)
                soln = solver.get_p()
                if prev_soln is not None:
                    assert_array_almost_equal(prev_soln, soln)
                prev_soln = soln
    
    test_support.run_unittest(ModelTest)
                


if __name__ == '__main__':
    test_models()