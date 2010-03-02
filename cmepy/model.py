"""
Model definition constants and validation functions
"""

PROPENSITIES = 'propensities'
TRANSITIONS = 'transitions'
NAME = 'name'
INITIAL_STATE = 'initial_state'
SHAPE = 'shape'
SPECIES_NAMES = 'species'
SPECIES_COUNTS = 'species_counts'
REACTION_NAMES = 'reactions'
ENTRIES = frozenset([
    PROPENSITIES,
    TRANSITIONS,
    NAME,
    INITIAL_STATE,
    SHAPE,
    SPECIES_NAMES,
    SPECIES_COUNTS,
    REACTION_NAMES,
])

class Model(dict):
    """
    Simple CME Model structure. Data can be accessed via key or attribute.
    """
    def __init__(self, **entries):
        dict.__init__(self, entries)
        self.validate()
    
    def __getattribute__(self, attrname):
        """
        allows keys to be accessed as attributes
        """
        try:
            return object.__getattribute__(self, attrname)
        except AttributeError:
            if attrname in self:
                return self[attrname]
            else:
                lament = 'model has no attribute \'%s\%'
                raise AttributeError(lament % str(attrname))
    
    def validate(self):
        """
        raise exception if model is not valid
        """
        validate_model(self)

def create(**entries):
    """
    Returns a model instance, using the specified entries.
    
    Entries should be model parameters as keyword arguments.
    """
    
    m = Model(**entries)
    m.validate()
    return m

def raise_error(error_type, descr, value, message):
    """
    raise_error(error_type, descr, value, message) -> raise error_type(...)
    """
    if descr is None:
        descr = 'object'
    lament = '%s of \'%s\' %s'
    raise error_type(
        lament % (
            descr,
            type(value),
            message
        )
    )

def must_be_within_range(value, descr=None, min_value=None, max_value=None):
    """
    raises ValueError if value is not within specified range
    """
    if (min_value is not None) and (value < min_value):
        message = 'has value %d < min %d' % (value, min_value)
        raise_error(ValueError, descr, value, message)
    if (max_value is not None) and (value > max_value):
        message = 'has value %d > max %d' % (value, max_value)
        raise_error(ValueError, descr, value, message)

def must_have_length(value, descr=None, min_len=None, max_len=None):
    """
    raises TypeError if value has no len()
    
    optionally checks to see if len(value) is inside specified range
    """
    if not hasattr(value, '__len__'):
        raise_error(TypeError, descr, value, 'has no len()')
    must_be_within_range(
        len(value),
        'len() of %s' % descr,
        min_len,
        max_len
    )

def must_be_callable(value, descr=None):
    """
    raises TypeError if value is not callable
    """
    if not hasattr(value, '__call__'):
        raise_error(TypeError, descr, value, 'is not callable')

def must_be_string_like(value, descr=None):
    """
    raises TypeError if value isn't string-like
    """
    try:
        str_value = str(value)
        if not value == str_value:
            raise ValueError
    except ValueError:
        raise_error(TypeError, descr, value, 'is not string-like')

def must_be_int_like(value, descr=None, min_value=None, max_value=None):
    """
    raises TypeError if value isn't int-like
    """
    try:
        str_value = int(value)
        if not value == str_value:
            raise ValueError
    except ValueError:
        raise_error(TypeError, descr, value, 'is not int-like')
    must_be_within_range(
        value,
        descr,
        min_value,
        max_value
    )

def must_be_mapping(value, descr=None):
    """
    raises TypeError if value isn't a mapping
    
    this is a very loose definition of a mapping.
    """
    if not hasattr(value, '__getitem__'):
        raise_error(TypeError, descr, value, 'is not a mapping')
    if not hasattr(value, '__contains__'):
        raise_error(TypeError, descr, value, 'is not a mapping')

def must_contain_item(value, descr, item):
    """
    raises KeyError if value doesn't contain item
    """
    if item not in value:
        message = 'does not contain: \'%s\'' % str(item)
        raise_error(KeyError, descr, value, message)

def validate_name(value):
    """
    raises exception if value is not a valid model name entry
    """
    must_be_string_like(value, 'name')

def validate_propensities(value):
    """
    raises exception if value is not a valid model propensities entry
    """
    must_have_length(value, 'propensity sequence', min_len=1)
    for prop in value:
        must_be_callable(prop, 'propensity function')

def validate_species_names(value):
    """
    raises exception if value is not a valid model species entry
    """
    must_have_length(value, 'species name sequence', min_len=1)
    for species_count in value:
        must_be_string_like(species_count, 'species name')

def validate_species_counts(value):
    """
    raises exception if value is not a valid model species_counts entry
    """
    must_have_length(value, 'species count sequence', min_len=1)
    for species_count in value:
        must_be_callable(species_count, 'species count function')

def validate_transitions(value):
    """
    raises exception if value is not a valid model transitions entry
    """
    must_have_length(value, 'transition sequence')
    for trans in value:
        must_have_length(trans, 'transition vector', min_len=1)
        for coord in trans:
            must_be_int_like(coord, 'transition vector element')

def validate_initial_state(value):
    """
    raises exception if value is not a valid model initial_state entry
    """
    must_have_length(value, 'initial state')
    for coord in value:
        must_be_int_like(coord, 'initial state element')

def validate_shape(value):
    """
    raises exception if value is not a valid model shape entry
    """
    must_have_length(value, 'shape')
    for coord in value:
        must_be_int_like(coord, 'shape element', min_value=1)
            
def validate_model(m):
    """
    raises exceptions if there's anything suspect about the model m
    
    this is not intended to be exhaustive or foolproof, it is intended merely
    to raise exceptions with more detailed and explanatory messages for the
    more obvious kinds of mistakes.
    """
    
    must_be_mapping(m, 'model')
    must_contain_item(m, 'model', PROPENSITIES)
    validate_propensities(m[PROPENSITIES])
    must_contain_item(m, 'model', TRANSITIONS)
    validate_transitions(m[TRANSITIONS])
    
    if len(m[PROPENSITIES]) != len(m[TRANSITIONS]):
        message = 'mismatched lengths of %s and %s in model %s'
        raise ValueError(
            message % (PROPENSITIES, TRANSITIONS, str(m))
        )
    
    if NAME in m:
        validate_name(m[NAME])
    if SPECIES_NAMES in m:
        validate_species_names(m[SPECIES_NAMES])
    if SPECIES_COUNTS in m:
        validate_species_counts(m[SPECIES_COUNTS])
    if INITIAL_STATE in m:
        validate_initial_state(m[INITIAL_STATE])
    if SHAPE in m:
        validate_shape(m[SHAPE])
    if REACTION_NAMES in m:
        validate_species_names(m[REACTION_NAMES])
    
    if (SPECIES_NAMES in m) and (SPECIES_COUNTS in m):
        if len(m[SPECIES_NAMES]) != len(m[SPECIES_COUNTS]):
            message = 'mismatched lengths of %s and %s in model %s'
            raise ValueError(
                message % (SPECIES_NAMES, SPECIES_COUNTS, str(m))
            )
    
    if (INITIAL_STATE in m) and (SHAPE in m):
        if len(m[INITIAL_STATE]) != len(m[SHAPE]):
            message = 'mismatched lengths of %s and %s in model %s'
            raise ValueError(
                message % (INITIAL_STATE, SHAPE, str(m))
            )
    if REACTION_NAMES in m:
        if len(m[REACTION_NAMES]) != len(m[PROPENSITIES]):
            message = 'mismatched lengths of %s and %s in model %s'
            raise ValueError(
                message % (REACTION_NAMES, PROPENSITIES, str(m))
            )
    
    for entry in m:
        if entry not in ENTRIES:
            message = 'has unexpected key \'%s\'' % entry
            raise_error(KeyError, 'model', m, message)
    
    # if we've got this far, model is probably hopefully okay-ish!
    return
