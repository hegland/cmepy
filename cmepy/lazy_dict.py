"""
Dictionary with lazy evaluation on access, via a supplied update function
"""

import itertools

class LazyDict(dict):
    """
    Dictionary that lazily updates values, when accessed, via supplied function.
    
    All the usual mapping methods should work, lazily updating, in a consistent
    fashion. Hopefully
    """
    def __init__(self, update_value, items = None):
        """
        LazyDict(update_value) -> lazy_dict
        
        Where
        
            update_value(key, existing_value, member) -> updated_value
        
        Suppose:
        
            lazy_dict[key] = value
        
        Then:
            v = lazy_dict[key]
            
                <=>
            
            lazy_dict[key] = update_value(key, value, (key in lazy_dict))
            v = update_value(key, value, (key in lazy_dict))
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
