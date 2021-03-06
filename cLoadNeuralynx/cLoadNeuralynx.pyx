'''
Wrapper of loadNeuralynx.c for loading Neuralynx data files.
'''

'''
Compile with:
python setup.py build_ext --inplace

Various resources can be found at:
http://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api

TO DO:
- How should I free the memory allocated to DataNCS.samples and DataNTT.samples?

'''

__version__ = '0.1'
__author__ = 'Santiago Jaramillo'
__created__ = '2011-01-13'

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free

np.import_array()  # Import the C-API for numpy

# -- For NCS files (continuous data) ---
# -- For NTT files (tetrode data) ---
cdef extern from "loadNeuralynx.h":
    ctypedef unsigned long long TStype
    ctypedef unsigned short SAMPLEtype
    ctypedef struct DataNCS:
        char* header
        TStype* timestamps
        unsigned long nRecords
        SAMPLEtype* samples
        unsigned long nSamples
    DataNCS readNCS(char* fileName)
    ctypedef struct DataNTT:
        char* header
        TStype* timestamps
        unsigned long nEvents
        SAMPLEtype* samples
        unsigned long nSamples
    DataNTT readNTT(char* fileName)


def read_ncs(fileName,readTS=False):
    '''Read Neuralynx NCS file containing continuous data.'''
    cdef DataNCS c_data     # Define the C structure to contain the data
    c_data = readNCS(fileName)
    nRecords = c_data.nRecords
    cdef np.npy_intp c_nRecords = <np.npy_intp>nRecords
    nSamples = c_data.nSamples
    cdef np.npy_intp c_nSamples = <np.npy_intp>nSamples
    samples = np.PyArray_SimpleNewFromData(1, &c_nSamples, np.NPY_SHORT,
                                           <void *>c_data.samples)
    if readTS:
        timestamps = np.PyArray_SimpleNewFromData(1, &c_nRecords, np.NPY_ULONGLONG,
                                                  <void *>c_data.timestamps)
        return (samples,timestamps)
    else:
        return samples


class DataCont(object):
    '''Access Neuralynx NCS files containing continuous data.

    FIXME: I'm currently making copies of the results of PyArray_SimpleNewFromData,
           because otherwise I don't know how to deallocate the memory of these arrays
           and I get a memory leak. Maybe the solution is related to the way I am
           allocating space on readNCS (loadNeuralynx.c).
    '''
    def __init__(self,fileName,readTS=True):
        cdef DataNCS c_data        # Define C structure to contain the data
        c_data = readNCS(fileName) # Read data from file
        self.filename = fileName
        self.header = c_data.header
        self.nRecords = c_data.nRecords
        self.nSamples = c_data.nSamples
        cdef np.npy_intp c_nSamples = <np.npy_intp>self.nSamples
        cdef np.npy_intp c_nRecords = <np.npy_intp>self.nRecords
        if readTS:
            self.timestamps = np.copy(np.PyArray_SimpleNewFromData(1, &c_nRecords, np.NPY_ULONGLONG,
                                                                   <void *>c_data.timestamps))
        else:
            self.timestamps = None
        self.samples = np.copy(np.PyArray_SimpleNewFromData(1, &c_nSamples, np.NPY_SHORT,
                                                            <void *>c_data.samples))
        free(c_data.header) # Deallocates, after a copy has been made
        free(c_data.timestamps) # Deallocates, after a copy has been made
        free(c_data.samples) # Deallocates, after a copy has been made

class DataTetrode(object):
    '''Access to Neuralynx NTT files containing tetrode data.

    FIXME: I'm currently making copies of the results of PyArray_SimpleNewFromData,
           because otherwise I don't know how to deallocate the memory of these arrays
           and I get a memory leak. Maybe the solution is related to the way I am
           allocating space on readNTT (loadNeuralynx.c).
    '''
    def __init__(self,fileName,readWaves=False):
        cdef DataNTT c_data        # Define C structure to contain the data
        c_data = readNTT(fileName) # Read data from file
        self.filename = fileName
        self.header = c_data.header
        self.nEvents = c_data.nEvents
        self.nSamples = c_data.nSamples
        cdef np.npy_intp c_nSamples = <np.npy_intp>self.nSamples
        cdef np.npy_intp c_nEvents = <np.npy_intp>self.nEvents
        self.timestamps = np.copy(np.PyArray_SimpleNewFromData(1, &c_nEvents, np.NPY_ULONGLONG,
                                                               <void *>c_data.timestamps))
        if readWaves:
            self.samples = np.copy(np.PyArray_SimpleNewFromData(1, &c_nSamples, np.NPY_SHORT,
                                                                <void *>c_data.samples))
        else:
            self.samples = None
        free(c_data.header) # Deallocates, after a copy has been made
        free(c_data.timestamps) # Deallocates, after a copy has been made
        free(c_data.samples) # Deallocates, after a copy has been made
        #self.timestamps = np.PyArray_SimpleNewFromData(1, &c_nEvents, np.NPY_ULONGLONG,
        #                                               <void *>c_data.timestamps)
        #tmp_timestamps = np.PyArray_SimpleNewFromData(1, &c_nEvents, np.NPY_ULONGLONG,
        #                                               <void *>c_data.timestamps)
        #self.timestamps = []#np.copy(tmp_timestamps)
        #del(tmp_timestamps) # Deletion of local or C global name not supported
        #free(c_data.timestamps) # Deallocates, but then there is no data (and segfault)
        #free(c_data) # Cannot assign type 'DataNTT' to 'void *'
        #np.PyArray_free(c_data.timestamps)  # DOES NOT WORK (type of pointer)

