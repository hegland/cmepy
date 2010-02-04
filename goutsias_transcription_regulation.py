"""
the transcription regulation example, as formulated by:

@article{goutsias2005quasiequilibrium,
  title={{Quasiequilibrium approximation of fast reaction
          kinetics in stochastic biochemical systems}},
  author={Goutsias, J.},
  journal={The Journal of chemical physics},
  volume={122},
  pages={184102},
  year={2005}
}

a time-independent approximation of this example was also considered by

@conference{burrage2006krylov,
  title={{A Krylov-based finite state projection algorithm for
          solving the chemical master equation arising in the
          discrete modelling of biological systems}},
  author={Burrage, K. and Hegland, M. and Macnamara, S. and Sidje, R.},
  booktitle={Proc. of The AA Markov 150th Anniversary Meeting},
  pages={21--37},
  year={2006}
}
"""

import numpy
import cmepy.solver
import cmepy.recorder
from cmepy.util import non_neg

def create_time_dependencies():
    """
    returns time dependencies as dictionary of functions keyed
    by sets of reaction indices
    """  
    
    avogadro_number = 6.0221415e23
    cell_v_0 = 1.0e-15 # litres
    k = avogadro_number*cell_v_0
    cell_t = 35*60.0 # seconds
    
    def phi(t):
        return 1.0/(k*numpy.exp(numpy.log(2.0)*t/cell_t))
    
    return {frozenset([4, 6, 8]) : phi}

def create_model(dna_count):

    c = {'DNA' : lambda *x : x[0],
         'DNA-D' : lambda *x : x[1],
         'DNA-2D' : lambda *x : dna_count - x[0] - x[1],
         'RNA' : lambda *x : x[2],
         'M' : lambda *x : x[3],
         'D' : lambda *x : x[4], }
    
    k = (4.3e-2,
         7.0e-4,
         7.8e-2,
         3.9e-3,
         1.2e7,
         4.791e-1,
         1.2e5,
         8.765e-12,
         1.0e8,
         0.5, )
    
    props = (lambda *x : k[0]*c['RNA'](*x),
             lambda *x : k[1]*c['M'](*x),
             lambda *x : k[2]*c['DNA-D'](*x),
             lambda *x : k[3]*c['RNA'](*x), 
             lambda *x : k[4]*c['D'](*x)*c['D'](*x), 
             lambda *x : k[5]*c['DNA-D'](*x), 
             lambda *x : k[6]*c['DNA-D'](*x)*c['D'](*x), 
             lambda *x : k[7]*c['DNA-2D'](*x), 
             lambda *x : k[8]*0.5*c['M'](*x)*non_neg(c['M'](*x)-1), 
             lambda *x : k[9]*c['D'](*x), )
    
    offsets = ((0, 0, 0, 1, 0),     # RNA -> RNA + M
               (0, 0, 0, -1, 0),    # M -> *
               (0, 0, 1, 0, 0),     # DNA-D -> RNA + DNA-D
               (0, 0, 1, 0, 0),     # RNA -> *
               (-1, 1, 0, 0, -1),   # DNA + D -> DNA-D
               (1, -1, 0, 0, 1),    # DNA-D -> DNA + D
               (0, -1, 0, 0, -1),   # DNA-D + D -> DNA-2D
               (0, 1, 0, 0, 1),     # DNA-2D -> DNA-D + D
               (0, 0, 0, -2, 1),    # M + M -> D
               (0, 0, 0, 2, -1), )  # D -> M + M
    
    species_names = []
    species_counts = []
    for key, value in c.iteritems():
        species_names.append(key)
        species_counts.append(value)
    
    origin = (dna_count, 0, 0, 2, 6)
    
    model = {'doc' : 'Goutsias transcription regulation',
             'species' : species_names,
             'species counts' : species_counts,
             'propensities' : props,
             'offset_vectors' : offsets,
             'norigin' : origin }
    return model

def main():
    dna_count = 2
    rna_max = 15
    m_max = 15
    d_max = 15
    model = create_model(dna_count)
    np = (dna_count+1, )*2 + (rna_max+1, m_max+1, d_max+1)
    model['np']= np
    
    solver = cmepy.solver.create(
        model,
        sink = True,
        time_dependencies = create_time_dependencies()
    )
    
    recorder = cmepy.recorder.create(
        ('species',
         model['species'],
         model['species counts'])
    )
    
    t_final = 0.01 
    time_steps = numpy.linspace(0.0, t_final, 21)
    for t in time_steps:
        print 't = %f' % t
        solver.step(t)
        p, p_sink = solver.y
        print 'p_sink = %f' % p_sink
        recorder.write(t, p)
    
    cmepy.recorder.display_plots(recorder, 'species', title = model['doc'])

if __name__ == '__main__':
    main()
