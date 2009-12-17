"""
CmeRecorder is a utility class for computation of measurements.
"""

import cmepy.new_core.sparse_marginal as sparse_marginal

def _clean_output_name(output_name):
    """
    used internally by MeasurementInfo
    """
    return output_name.replace(' ', '_')

class MeasurementInfo(object):
    """
    Stores measurements of a specified variable.
    
    Created internally by CmeRecorder.
    """
    
    def __init__(self, name, target_outputs):
        """
        MeasurementInfo instances are constructed internally via CmeRecorder
        
        Arguments:
        
        name : name of the variable
        target_outputs : list of output types to store
        """
        
        self.name = name
        self.target_outputs = target_outputs
        for output in self.target_outputs:
            self.__dict__[_clean_output_name(output)] = []
        
        self._transform = None
        self._dim = None
        self.times = []
    
    def set_transform(self, transform):
        """
        Specifies how to transform from the state space of the probability
        distribution p to the state space of this variable.
        
        Mutually exclusive to set_dim.
        """
        if self._dim is not None:
            raise RuntimeError('cannot set both transform and dim')
        self._transform = transform
    
    def set_dim(self, dim):
        """
        Specifies which dimension of the state space this variable corresponds
        to.
        
        Mutually exclusive to set_transform.
        """
        if self._transform is not None:
            raise RuntimeError('cannot set both transform and dim')
        self._dim = dim
    
    def compute_marginal(self, p, norigin):
        if self._transform is not None:
            marginal = sparse_marginal.create_sparse_marginal(self._transform,
                                                              p,
                                                              norigin)
        elif self._dim is not None:
            marginal = sparse_marginal.create_coord_sparse_marginal(self._dim,
                                                                    p,
                                                                    norigin)
        else:
            raise RuntimeError('transform or dim must first be set')
        return marginal
    
    def add_measurement(self, output_name, measurement):
        """
        Stores given measurement of specified output.
        
        Called via CmeRecorder to add measurement entries.
        """
        self.__dict__[_clean_output_name(output_name)].append(measurement)

class CmeRecorder(object):
    """
    CmeRecorder is a utility class to compute common measurements,
    such as marginals, expected values, and standard deviations, from
    a given disitribution p. 
    """
    
    def __init__(self, model, **kwargs):
        """
        Arguments:
        
        model : a CME model
        
        Optional Keyword Arguments:
        
        norigin : origin of the state space used.
            Defaults to (0,0, ..., 0).
        """
        self.model = model
        self.args = kwargs
        if 'norigin' not in self.args:
            self.args['norigin'] = None
        
        self._measurements = {}
        self.group_names = set()
        self.output_names = set(['marginal',
                                 'expected value',
                                 'standard deviation'])
        self._measurement_cache = None
        
    def get_measurement_info(self, group_name, var_name):
        """
        Returns recorded measurement information for specified variable from
        specified group.
        """
        return self._measurements[group_name][var_name]
    
    def add_target(self, group_name, outputs, variable_names, transforms=None):
        """
        Register output types to be computed for a group of variables.
        """
        
        for output_name in outputs:
            if output_name not in self.output_names:
                complaint = ('unknown output format: \"'+str(output_name)
                             +'\". Known formats: ')
                known_output_formats = ', '.join(self.output_names)
                raise KeyError, complaint + known_output_formats
        
        if group_name not in self.model:
            complaint = ('model contains no entry for the group name'
                         +' \"'+str(group_name)+'\".')
            raise KeyError, complaint
        
        self.group_names.add(group_name)
        
        if transforms is not None:
            if len(transforms) != len(variable_names):
                complaint = 'lengths of transforms and variable_names disagree'
                raise ValueError, complaint
        
        for i, name in enumerate(variable_names):
            # compute dim, the dimension this name corresponds to
            # we do this by linearly searching through the names listed
            # by the model for this group (reaction / species).
            model_var_names = self.model[group_name]
            var_dim = None
            for dim in xrange(len(model_var_names)):
                if model_var_names[dim] == name:
                    var_dim = dim
                    break
            if var_dim is None:
                complaint = (('model[\'%s\'] does not contain'%group_name) +
                             (' the variable name \'%s\''%name))
                raise KeyError, complaint 
            
            m_info = MeasurementInfo(name, outputs)
            
            if transforms is not None:
                m_info.set_transform(transforms[i])
            else:
                m_info.set_dim(var_dim)
            
            if group_name not in self._measurements:
                self._measurements[group_name] = {}
            self._measurements[group_name][name] = m_info
                                                                  
                                                                  
    
    def _compute(self,
                group_name,
                output_name,
                var_name,
                p_current,
                t_current):
        
        key = (group_name, output_name, var_name)
        m_info = self.get_measurement_info(group_name, var_name)
        
        # if we've previously computed this measurement, just return that!
        if key in self._measurement_cache:
            return self._measurement_cache[key]
        else:
            # otherwise, we need to compute this measurement
            
            # determine the prerequisite measurements required to
            # compute this one
            marginal_key = (group_name,
                            'marginal',
                            var_name)
            expected_value_key = (group_name,
                                  'expected value',
                                  var_name)
            
            if output_name == 'marginal':
                prereq_keys = set()
            if output_name == 'expected value':
                prereq_keys = set([marginal_key])
            elif output_name == 'standard deviation':
                prereq_keys = set([marginal_key,
                                   expected_value_key])
            
            
            # check to see if the prerequisites were already computed - if they
            # weren't, we recursively compute them
            for prereq_key in prereq_keys:
                if prereq_key not in self._measurement_cache:
                    prereq_group, prereq_output, prereq_var = prereq_key
                    self._compute(prereq_group,
                                  prereq_output,
                                  prereq_var,
                                  p_current,
                                  t_current)
                
                assert (prereq_key in self._measurement_cache)
            
            # since all the pre-reqs are computed, we can now compute this
            # measurement
            if output_name == 'marginal':
                marginal = m_info.compute_marginal(p_current,
                                                   self.args['norigin'])
                measurement = marginal
            elif output_name == 'expected value':
                marginal = self._measurement_cache[marginal_key]
                expected_value = marginal.expected_value()
                measurement = expected_value
            elif output_name == 'standard deviation':
                marginal = self._measurement_cache[marginal_key]
                expected_value = self._measurement_cache[expected_value_key]
                standard_deviation = marginal.standard_deviation(expected_value)
                measurement = expected_value = standard_deviation
        
        self._measurement_cache[key] = measurement
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
        self._measurement_cache = {}
        for group_name in self.group_names:
            for var_name in self._measurements[group_name]:
                m_info = self._measurements[group_name][var_name]
                m_info.times.append(t)
                for output_name in m_info.target_outputs:
                    measurement = self._compute(group_name,
                                                output_name,
                                                var_name,
                                                p,
                                                t)
                    m_info.add_measurement(output_name, measurement)
        self._measurement_cache = None
                
    
    def measurements(self, group_name):
        """
        Generator yielding the measurement information for all variables
        from the specified group.
        
        Arguments:
        
        group_name : name of a group previously mentioned by a add_target
            call.
        """
        for var_name in self._measurements[group_name]:
            measurement_info = self._measurements[group_name][var_name]
            yield measurement_info
