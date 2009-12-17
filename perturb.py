

def get_simple_model(count=1):
    epsilon = 0.01
    epsilon_inverse = 1.0 / epsilon
    model = {'propensities' : (lambda *x : 1.0*x[0],
                               lambda *x : epsilon_inverse*x[1]),
             'np' : (count+1, )*2,
             'offset_vectors' : ((-1, 1), (0, -1)),
             'doc' : 'simple catalytic reaction scheme',
             'species' : ('S_1', 'S_2', 'S_3'),
             'species counts' : (lambda *x : x[0],
                                 lambda *x : x[1],
                                 lambda *x : count-x[0]-x[1])
             }
    return model

def compute_system_matrix(model):
    from cmepy.core.matrix_cme import (process_offset_vectors,
                                       gen_reaction_actions,
                                       gen_sparse_matrix)
    np = model['np']
    propensities = model['propensities']
    norigin = model.get('norigin', None)
    offset_vectors = model.get('offset_vectors', None)
    offset_vectors = process_offset_vectors(np,
                                            propensities,
                                            offset_vectors) 
    reaction_actions = gen_reaction_actions(np,
                                            propensities,
                                            offset_vectors,
                                            norigin)
    matrix = gen_sparse_matrix(np, reaction_actions)
    return matrix

def decompose_model(full_model, reaction_subsets):
    from itertools import izip
    
    def partially_copy_model(model):
        model_copy = dict(model)
        keys_to_remove = ('propensities',
                          'offset_vectors',
                          'reaction_names', )
        for key in keys_to_remove:
            if key in model_copy:
                del model_copy[key]
        
        return model_copy
    
    def append_entry(model, key, entry):
        entries = model.get(key, [])
        entries.append(entry)
        model[key] = entries
    
    subsets_to_models = {}
    for subset in reaction_subsets:
        subsets_to_models[subset] = partially_copy_model(full_model)
    
    for i, (prop, offset) in enumerate(izip(full_model['propensities'],
                                            full_model['offset_vectors'])):
        for subset in subsets_to_models:
            if i in subset:
                subset_model = subsets_to_models[subset]
                append_entry(subset_model, 'propensities', prop)
                append_entry(subset_model, 'offset_vectors', offset)
    
    return tuple(subsets_to_models.values())

def truncated_svd(a, k):
    """
    approximate svd involving only the k largest singular values
    """
    import scipy.linalg
    
    u, s, vh = scipy.linalg.svd(a)
    u_bar = u[:, :k]
    s_bar = s[:k]
    vh_bar = vh[:k, :]
    
    return (u_bar, s_bar, vh_bar)
    