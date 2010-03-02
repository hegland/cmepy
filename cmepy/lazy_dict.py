"""
Dictionary with lazy evaluation on access, via a supplied update function
"""

import itertools

class LazyDict(dict):
    """
    A dictionary type that lazily updates values when they are accessed.
    
    All the usual dictionary methods work as expected, with automatic lazy
    updates occuring behind the scenes whenever values are read from the
    dictionary. 
    
    The optional ``items`` argument, if specified, is a mapping instance used
    to initialise the items in the :class:`LazyDict`.
    
    The ``update_value`` argument required by the :class:`LazyDict` constructor
    must be a function of the form:
    
        update_value(k, existing_value, member) -> updated_value
    
    This function is called whenever an item with the key ``k`` is read
    from the :class:`LazyDict`. The second argument ``existing_value``, is
    the value corresponding to the key ``k`` stored in the :class:`LazyDict`,
    or ``None``, if the key ``k`` is not contained in the :class:`LazyDict`.
    The third argument ``member`` is a boolean value indicating if there is
    an existing value stored under the key ``k``.
    
    This function is used as follows by the :class:`LazyDict`. Suppose that the
    value ``v`` has been stored in a :class:`LazyDict` object ``lazy_dict``
    under the key ``k``, that is, ``lazy_dict[k] = v``. Then subsequently
    accessing this value in the usual manner::
    
        v_updated = lazy_dict[k]
    
    is equivalent to the following two statements::
    
        lazy_dict[k] = update_value(k, v, (k in lazy_dict))
        v_updated = update_value(k, v, (k in lazy_dict))
    
    Observe how the value stored in the :class:`LazyDict` under the key ``k``
    is first updated, using the provided function,
    with the updated value then being the one returned.
    """
    def __init__(self, update_value, items = None):
        """
        Returns a LazyDict using the specified ``update_value`` function
        and optional initial dictionary arguments.
        """
        self.update_value = update_value
        if items is None:
            dict.__init__(self)
        else:
            dict.__init__(items)
        
    def __getitem__(self, key):
        member = dict.__contains__(self, key)
        if member:
            existing_value = dict.__getitem__(self, key)
        else:
            existing_value = None
        # ensure measurement is up to date
        updated_value = self.update_value(key, existing_value, member)
        self[key] = updated_value
        return updated_value
    
    def copy(self):
        return LazyDict(self.update_value, dict.copy(self))
    
    def itervalues(self):
        return itertools.imap((lambda k : self[k]), dict.iterkeys(self))
    
    def iteritems(self):
        return itertools.imap((lambda k : (k, self[k])), dict.iterkeys(self))
    
    def pop(self, *args):
        n_args = len(args)
        if n_args < 1:
            raise TypeError('pop expected at least 1 argument, got %d' % n_args)
        if n_args > 2:
            raise TypeError('pop expected at most 2 arguments, got %d' % n_args)
        
        k = args[0]    
        if k in self:
            value = self[k]
            del self[k]
            return value
        else:
            if n_args == 2:
                return args[1]
            else:
                raise KeyError(str(k))
    
    def popitem(self):
        key, value = dict.popitem(self)
        self[key] = value
        updated_value = self[key]
        del self[key]
        return key, updated_value
    
    def setdefault(self, k, x=None):
        if k in self:
            return self[k]
        else:
            self[k] = x
            return x
    
    def get(self, k, x=None):
        if k in self:
            return self[k]
        else:
            return x
    
    def values(self):
        return list(self.itervalues())
    
    def items(self):
        return list(self.iteritems())
