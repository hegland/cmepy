"""
CmeRecorder is a utility class for computation of measurements.
"""

import cmepy.new_core.sparse_marginal as sparse_marginal


class MeasurementInfo(object):
    """
    Stores measurements of a specified variable.
    
    Created internally by CmeRecorder.
    """
    
    def __init__(self, name, dim, target_outputs):
        """
        MeasurementInfo instances are constructed internally via CmeRecorder
        
        Arguments:
        
        name : name of the variable
        dim : dimension index corresponding to this name (TODO FIXME ...)
        target_outputs : sequence of outputs to store for this variable
        """
        self.name = name
        self.dim = dim
        self.times = []
        self.target_outputs = target_outputs
        for output in target_outputs:
            self.__dict__[output] = []
    
    def add_measurement(self, output_name, measurement):
        """
        Stores given measurement of specified output.
        
        Called via CmeRecorder to add measurement entries.
        """
        self.__dict__[output_name].append(measurement)

class Measurer(object):
    
    def __init__(self, group_name, variable_name, marginal_creators):
        self.group_name = group_name
        self.output_name = output_name
        self.variable_names = variable_names
        self.marginal_creators = marginal_creators
    
    def compute(self, var_name, t, p, measurement_cache):

class CmeRecorder(object):
    """
    XXX TODO
        
    """
    
    def __init__(self, model, **kwargs):
        """
        Arguments:
        
        model : a CME model
        
        Optional Keyword Arguments:
        
        norigin : origin of the reaction-count state space used by the
            CmeSolver. Defaults to (0,0, ..., 0).
        """
        self.model = model
        self.dims = len(self.model['np'])
        self.args = kwargs
        if 'norigin' not in self.args:
            self.args['norigin'] = (0,)*self.dims
        
        self.measurers = {}
        self.group_names = set()
        self.output_names = set()
        self.__measurement_cache = None
    
    def __add_measurer(self, group_name, output_name, variable_name, transform):
        
        
    
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
                if group_name not in self.model:
                    complaint = ('model contains no variable names for the key'
                                 +' \"'+str(group_name)+'\".')
                    raise KeyError, complaint
                model_var_names = self.model[group_name]
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
                
    
    def write(self, t, p):
        """
        Computes the measurements, specified previously via add_target,
        using the given cme probability distribution p and time t.
        """
        
        # accumulate a cache of measurements
        # some measurements have pre-req measurements, which in general
        # are not wanted by the user, so we cache them here for the
        # duration of this call
        self.__measurement_cache = {}
        for group_name in self.group_names:
            for var_name in self.__dict__[group_name]:
                m_info = self.__dict__[group_name][var_name]
                m_info.times.append(t)
                for output_name in m_info.target_output:
                    measurement = self.__compute(group_name,
                                                 output_name,
                                                 var_name,
                                                 p,
                                                 t)
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
