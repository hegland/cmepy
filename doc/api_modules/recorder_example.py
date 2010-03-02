"""
example of basic recorder usage
"""


def main():
    """
    plot example graph of species count expected value over time
    """
    from cmepy import recorder
    
    # create recorder, adding target random variable
    # for species counts of A and B
    r = recorder.create(
        (['A', 'B'], [lambda *x : x[0], lambda *x : x[1]])
    )
    
    # write some solution data to the recorder
    r.write(0, {(1, 0) : 1.0})
    r.write(0.5, {(1, 0) : 0.9, (0, 1) : 0.1})
    
    # get measurements for species counts of B
    measurement = r['B']
    
    import pylab
    
    # plot expected value of species counts of B over time
    pylab.plot(measurement.times, measurement.expected_value, 'o')
    
    
    pylab.xlim(-0.1, 0.6)
    pylab.xlabel('solution time $t$')
    pylab.ylim(-0.1, 1.1)
    pylab.ylabel('expected value of species count $B$')
    pylab.title('expected value of species count $B$ over time')
    pylab.show()

if __name__ == '__main__':
    main()
