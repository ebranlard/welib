""" 
Set of functions to manipulate dictionaries
"""

class DictLibError(Exception):
    pass

def renameDict(d, keyMap, full_map=False, keys_exist=False):
    """ 
    Rename keys of a dictionary, preserving "order" of keys

    Missing keys

    INPUTS:
     - d: input dictionary
     - keyMap: dictionary mapping old keys to new keys for instance: {'old1':'new1', 'old2':'new2'}
     - full_map: if True, keyMap and d have exactly the same keys (possibly in different order).
     - keys_exist: if True, the keys of keyMap must be keys that are present in d (a subset or full).
    """

    # Sanity checks if user require some strict renaming
    S1 =  set(d.keys())
    S2 =  set(keyMap.keys())
    diff1 = S1.difference(S2) # Present in 1, absent in 2
    diff2 = S2.difference(S1) # Present in 2, absent in 1
    if full_map and S1 != S2:
        raise DictLibError('Keys of input dictionary are not all in keyMap but full_map is True.\nKey present in input but absent in map: {}.\nKeys present in map but absent in input: {}'.format(diff1, diff2))
    if keys_exist and len(diff2)>0:
        raise DictLibError('Keys of keyMap dictionary are not all in input dictionary but `keys_exist` is True.\nKeys present in map but absent in input: {}'.format(diff1, diff2))

    # Create new dictionary (preserves "order")
    d_out = {}
    for k,v in d.items():
        if k in keyMap.keys():
            k_new = keyMap[k]
            d_out[k_new] = v
        else:
           d_out[k] = v
    return d_out


def renameDictKey(d, old_key_pattern, new_key, regexp=False, key_exist=False, case_sensitive=False):
    """ 
    Change the key of a dictionary to `new_key`, if the old key is anything like possible_old_keys

    INPUTS:
     - d: input dictionary
     - old_key_pattern: either:
          - old key (string)
          - old key regexp pattern (string)
          - list of possible old keys, to be replaced by new_key
          - list of possible old regexp patterns, to be replaced by new_key
     - regexp: if True, old_key_pattern is a regexp 
     - key_exist: if True, one of the old_keys must be present in d.

    """
    # Sanity of inputs
    if isinstance(old_key_pattern, str):
        old_key_pattern=[old_key_pattern]
    if regexp:
        raise NotImplementedError() # ...
    if not case_sensitive:
        old_key_pattern=list(set([k.lower() for k in old_key_pattern]))

    # Create new dictionary (preserves "order")
    d_out = {}
    replaced = False
    for k,v in d.items():
        if not case_sensitive:
            k_old = k.lower()
        else:
            k_old = k

        if k_old in old_key_pattern:
            if replaced is True:
                raise DictLibError('Multiple keys of input dictionary match. This is not allowed.\nInput keys:{}\nOld keys pattern:{}'.format(d.keys(), old_key_pattern))
            d_out[new_key] = v
            replaced=True
        else:
            d_out[k] = v

    if not replaced and key_exist:
        raise DictLibError('Unable to find a match in dictionary for {}\n. The keys present are: {}'.format(old_key_pattern, d.keys()))

    return d_out






