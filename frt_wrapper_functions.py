#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 12:15:32 2019
Some helper functions to handle data moving between fortran and python.
Used by FRT to read input from input files which are not json.
@author: Matti Möttus
"""
import numpy as np

def arraychars( array_in, bufferlength=None ):
    """ Convert the array of strings in array_in into an array of characters
    required by F2PY. Adds one dimension to the array.
    bufferlength: length of strings as defined in Fortran"""

    if bufferlength is None:
        # attempt to recover it from the data type of array_in
        #    -- assuming it is an array of fixed-length string
        bufferlength = array_in.dtype.itemsize

    array_out = np.empty( list(array_in.shape) + [bufferlength], dtype='c' )
    # convert recursively, starting from the outmost dimension
    arraychars_recursive(array_in,array_out,bufferlength)
    return array_out

def arraychars_recursive( array_in, array_out, N ):
    """ Recursive helper function for arraychars
        array_in is input, array_out output
        N is maximum string length
    """
    # check array_in dimensionality to see if we need to recurse
    #  note: array_out has one more dimension compared with array_in
    if array_in.ndim > 0:
        for a_in,a_out in zip( array_in, array_out ):
            arraychars_recursive(a_in,a_out,N)
    else:
        # array_in is actually a scalar, i.e. a string
        # we are at the last index, fill the array element
        array_out[:] = chr(0) # set to blank
        str_out = array_in.astype(str)
        NN = min((N,len(str_out)))
        array_out[0:NN] = str_out[0:NN]

def str2chararray( str_in, N ):
    """ Convert a string into array of characters so they can be assigned
    to strings of fixed length N in Fortran common blocks using f2py """
    array_out = np.empty( N, dtype='c' )
    array_out[:] = chr(0) # set to blank
    NN = min((N,len(str_in)))
    array_out[0:NN] = str_in[0:NN]
    return array_out

def chararray2strarray( S ):
    """ convert the arrays of characters returned by fortran functions into arrays of string
    which can be easily manipulated in python.
    More explicitly: if a fortran function returns a "character*N S(M)", f2py retreives
    it as an M-element array of N-element array of characters. This function converts it to an
    M-element array of strings.
    """
    M = len(S[0])
    S_2 = S.reshape(-1,M).view( 'S'+str(M) )
    return [ x[0].decode('utf8').strip() for x in S_2 ]

