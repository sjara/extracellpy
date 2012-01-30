#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Extra functions that do not exist yet in Numpy
'''

import numpy as np

def index_to_bool(arrayOfIndexes):
    '''Convert array of indexes into a 2D array,
    where each column indicates if element'''
    elements = np.unique(arrayOfIndexes)
    arrayOfBool = np.zeros((arrayOfIndexes.shape[0],elements.shape[0]),dtype=bool)
    for inde,elem in enumerate(elements):
        arrayOfBool[:,inde] = arrayOfIndexes==elem
    return(arrayOfBool,elements)

def bool_to_index(arrayOfBool):
    pass

'''From http://projects.scipy.org/numpy/ticket/1036'''
def ismember_newer(totest, members):
    """
    A setmember1d, which works for totest arrays with duplicate values
    """
    uniques_in_test, rev_idx = np.unique1d(totest, return_inverse=True)
    uniques_in_members_mask = np.setmember1d(uniques_in_test, members)
    # Use this instead if members is not unique
    # uniques_in_members_mask = setmember1d(uniques_in_test, unique1d(members))
    return uniques_in_members_mask[rev_idx]

def ismember(a1,a2):
    """ Test whether items from a2 are in a1.
    
    This does the same thing as np.setmember1d, but works on
    non-unique arrays.
    
    Only a few (2-4) times slower than np.setmember1d, and a lot
    faster than [i in a2 for i in a1].

    An example that np.setmember1d gets wrong: 
    
    >>> a1 = np.array([5,4,5,3,4,4,3,4,3,5,2,1,5,5])
    >>> a2 = [2,3,4]
    >>> mask = ismember(a1,a2)
    >>> a1[mask]
    array([4, 3, 4, 4, 3, 4, 3, 2])
    """
    a2 = set(a2)
    a1 = np.asarray(a1)
    ind = a1.argsort()
    a1 = a1[ind]
    mask  = []
    # need this bit because prev is not defined for first item
    item  = a1[0]
    if item in a2:
        mask.append(True)
        a2.remove(item)
    else:
        mask.append(False)
    prev = item
    # main loop
    for item in a1[1:]:
        if item == prev:
            mask.append(mask[-1])
        elif item in a2:
            mask.append(True)
            prev = item
            a2.remove(item)
        else:
            mask.append(False)
            prev = item
    # restore mask to original ordering of a1 and return
    mask = np.array(mask)
    return mask[ind.argsort()]