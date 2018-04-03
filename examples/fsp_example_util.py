"""
common utility functions for the fsp example scripts
"""

def plot_solution_and_domain(measurement, domains):
    """
    displays plots of the solution and the domain at various times
    
    used by the three three fsp example scripts
    """
    import pylab
    # plot the solution
    shape = (41, 41)
    for i, (t, distribution) in enumerate(zip(measurement.times, measurement.distributions)):
        pylab.subplot(3, 3, i + 1)
        dense_distribution = distribution.to_dense(shape)
        pylab.imshow(
            dense_distribution,
            interpolation = 'nearest',
            origin = 'lower'
        )
    pylab.figure()
    # plot the states in the domain
    for i, (t, domain) in enumerate(zip(measurement.times, domains)):
        pylab.subplot(3, 3, i + 1)
        domain_x, domain_y = domain
        pylab.scatter(domain_x, domain_y, marker = 'o', c = 'k', s = 6)
        pylab.xlim(0, 45)
        pylab.ylim(0, 45)
    pylab.show()
