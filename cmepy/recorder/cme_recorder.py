"""
CmeRecorder is a utility class for computation of measurements.
"""

class CmeRecorder():
    """
    Computes and stores measurements from a CmeSolver instance.
    
    See documentation for the add_target method for a discussion
    of which variables & outputs can be measured.
    
    Example of usage:
        
        # initialise CmeSolver instance from some model
        solver = cmepy.solver.CmeSolver(model)
        
        # initialise CmeRecorder
        recorder = cmepy.recorder.CmeRecorder(solver)
        
        # specify which outputs we are interested in
        # e.g. here we wish to record expectation values
        # for all reactions named by the model
        recorder.add_target(output = ['expectation'],
                            reactions = model['reactions'])
        
        # advance solution of CME, taking measurements
        for t in time_steps:
            solver.step(t)
            recorder.take_measurements()
        
        # plot expected values of reaction counts
        import pylab
        pylab.figure()
        for r_info in recorder.measurements('reactions'):
            pylab.plot(r_info.times, # plot times measurements were recorded
                r_info.expectation,  # against expectation values
                label = r_info.name) # label graphs by reaction name
        pylab.legend()
        pylab.title('reaction count expected values')
        
    """
    
    class MeasurementInfo():
        """
        Stores measurements of a specified variable.
        
        Created internally by CmeRecorder.
        """
        
        def __init__(self, name, dim, target_output):
            """
            MeasurementInfo instances are constructed internally via CmeRecorder
            
            Arguments:
            
            name : name of the variable
            dim : dimension index corresponding to this name (TODO FIXME ...)
            target_output : sequence of outputs to store for this variable
            """
            self.name = name
            self.dim = dim
            self.times = []
            self.target_output = target_output
            for output in target_output:
                self.__dict__[output] = []
        
        def add_measurement(self, output_name, measurement):
            """
            Stores given measurement of specified output.
            
            Called via CmeRecorder to add measurement entries.
            """
            self.__dict__[output_name].append(measurement)
    
    def __add_all_measurers(self):
        """
        This is where we register all of the defined measurer subclasses ...
        """
        import cmepy.recorder.reaction_measurers as reaction_m
        import cmepy.recorder.species_measurers as species_m
        measurer_modules = [reaction_m, species_m]
        for measurer_module in measurer_modules:
            for measurer_class in measurer_module.all_measurers():
                measurer = measurer_class(self)
                self.__add_measurer(measurer)
    
    def __init__(self, cme_solver, **kwargs):
        """
        Arguments:
        
        cme_solver : cmepy.solver.CmeSolver instance from which to
                     take measurements.
        
        Optional Keyword Arguments:
        
        norigin : origin of the reaction-count state space used by the
            CmeSolver. Defaults to (0,0, ..., 0).
        """
        self.solver = cme_solver
        self.dims = len(self.solver.model['propensities'])
        self.args = kwargs
        if 'norigin' not in self.args:
            self.args['norigin'] = (0,)*self.dims
        
        self.measurers = {}
        self.group_names = set()
        self.output_names = set()
        self.__measurement_cache = None
        self.__add_all_measurers()
    
    def __add_measurer(self, measurer):
        key = (measurer.group_name, measurer.output_name)
        self.measurers[key] = measurer
        
        if measurer.group_name not in self.group_names:
            self.__dict__[measurer.group_name] = {}
            self.group_names.add(measurer.group_name)
        if measurer.output_name not in self.output_names:
            self.output_names.add(measurer.output_name)
    
    def get_measurement_info(self, group_name, var_name):
        """
        Returns recorded measurement information for specified variable from
        specified group.
        """
        return self.__dict__[group_name][var_name]
    
    def add_target(self, output, **kwargs):
        """
        Add specified outputs & variables for later measurement.
        
        Arguments:
        
        output : a sequence of valid output names.
            Valid output names for reactions:
                'expectation', 'std_dev', 'marginal'
            Valid output names for species:
                'expectation', 'std_dev', 'sparse_marginal'
        
        Optional Keyword Arguments:
        
        reactions : a sequence of reaction names defined by the model
        species : a sequence of species names defined by the model
        
        Some examples of usage:
        
            i.    To specify that the expectation of all reaction counts
                  should be measured:

                      recorder.add_target(output = ['expectation'],
                                          reactions = model['reactions'])
            
            ii.    To specify that the expectation and standard deviation
                   of all species counts should be measured:

                      recorder.add_target(output = ['expectation', 'std_dev'],
                                          species = model['species'])
            
            iii.   To specify that the marginal distribution of the reaction
                   count with name 'A->B'should be measured:

                      recorder.add_target(output = ['marginal'],
                                          reactions = ['A->B'])
        """
        
        for output_name in output:
            if output_name not in self.output_names:
                complaint = ('unknown output format: \"'+str(output_name)
                             +'\". Known formats: ')
                known_output_formats = ', '.join(self.output_names)
                raise KeyError, complaint + known_output_formats
        for group_name in kwargs:
            for name in kwargs.get(group_name, []):
                # compute dim, the dimension this name corresponds to
                # we do this by linearly searching through the names listed
                # by the model for this group (reaction / species).
                if group_name not in self.__dict__:
                    complaint = ('unknown variable group'
                                 +' \"'+str(group_name)+'\".'
                                 +' Known variable groups: ')
                    known_group_names = ', '.join(self.group_names)
                    raise KeyError, complaint + known_group_names
                if group_name not in self.solver.model:
                    complaint = ('model contains no variable names for the key'
                                 +' \"'+str(group_name)+'\".')
                    raise KeyError, complaint
                model_var_names = self.solver.model[group_name]
                var_dim = None
                for dim in xrange(len(model_var_names)):
                    if model_var_names[dim] == name:
                        var_dim = dim
                        break
                if var_dim is None:
                    raise Exception, 'unable to find name '+str(name)
                
                m_info = CmeRecorder.MeasurementInfo(name, var_dim, output)
                self.__dict__[group_name][name] = m_info
                                                                  
                                                                  
    
    def __compute(self,
                  group_name,
                  output_name,
                  var_name,
                  p_current,
                  t_current):
        
        key = (group_name, output_name, var_name)
        # if we've previously computed this measurement, just return that!
        if key in self.__measurement_cache:
            return self.__measurement_cache[key]
        else:
            # otherwise, we need to compute this measurement
            
            # look up the measurer instance which we'll use to
            # compute this measurements
            measurer = self.measurers[(group_name, output_name)]
            
            # determine the pre-requisite measurements required to
            # compute this one, by querying the measurer instance
            prereq_keys = measurer.get_prereq_keys(var_name)
            
            # check to see if the pre-reqs were already computed - if they
            # weren't, we recursively compute them
            for prereq_key in prereq_keys:
                if prereq_key not in self.__measurement_cache:
                    prereq_group, prereq_output, prereq_var = prereq_key
                    self.__compute(prereq_group,
                                   prereq_output,
                                   prereq_var,
                                   p_current,
                                   t_current)
                
                assert (prereq_key in self.__measurement_cache)
            
            # since all the pre-reqs are computed, we can now use the
            # measurer to compute this measurement
            # n.b. measurement_cache is passed, which contains all the
            # pre-req measurements required for this operation.
            measurement = measurer.compute(var_name,
                                           t_current,
                                           p_current,
                                           self.__measurement_cache)
            self.__measurement_cache[key] = measurement
        return measurement
                
    
    def take_measurements(self):
        """
        Computes the measurements, specified previously via add_target,
        using the current cme probability distribution p and time t
        obtained from the solver.
        """
        p_current = self.solver.get_p()
        t_current = self.solver.t
        
        # accumulate a cache of measurements
        # some measurements have pre-req measurements, which in general
        # are not wanted by the user, so we cache them here for the
        # duration of this call
        self.__measurement_cache = {}
        for group_name in self.group_names:
            for var_name in self.__dict__[group_name]:
                m_info = self.__dict__[group_name][var_name]
                m_info.times.append(t_current)
                for output_name in m_info.target_output:
                    measurement = self.__compute(group_name,
                                                 output_name,
                                                 var_name,
                                                 p_current,
                                                 t_current)
                    m_info.add_measurement(output_name, measurement)
        self.__measurement_cache = None
                
    
    def measurements(self, group_name):
        """
        Generator yielding the measurement information for all variables
        from the specified group.
        
        Arguments:
        
        group_name : name of a group previously mentioned by a add_target
            call.
        """
        for var_name in self.__dict__[group_name]:
            measurement_info = self.__dict__[group_name][var_name]
            yield measurement_info
