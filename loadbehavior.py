#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Read behavior data'''

import numpy as np

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

import tables
import os
import string
import struct

#BehavDataDir = '/var/data/test/'
#BehavFileName = os.path.join(BehavDataDir,'data_saja_reversal_santiago_saja064_20110204a.h5')

def pp(rawev):
    '''Pretty print a matrix of raw events
    It displays only four of the five columns, since the last one is irrelevant'''
    for row in rawev:
        print '%d \t % d \t %0.4f \t %d'%(row[0],row[1],row[2],row[3])

def ints_to_string(intarray):
    '''Convert array of integers to array of strings
    This method is currently very slow (because strings are stored for each trial)
    USE struct.pack INSTEAD.
    '''
    newarray = []
    for row in intarray:
        onestring=''.join([chr(c) for c in row])
        newarray.append(string.rstrip(onestring))
    return newarray


FORMAT_VERSION_HDF5_MATLAB = '1.3_BControl'
class BehaviorData(dict):
    '''Access to behavior data in HDF5 format saved by SaveToHDF5.m'''

    def __init__(self,behavFileName):
        h5file = tables.openFile(behavFileName, mode = "r")
        self.formatversion = h5file.getNode('/DataFormatVersion').read()[0]
        self.check_version()
        self.datanode = h5file.getNode('/SessionParams')
        self.populate_dict()
        h5file.close()

        # -- Change first index to zero (the Python/C way) --
        self['StateLabelTrialID'] = self['StateLabelTrialID'].astype('int32') - 1
        self['RawEventsTrialID'] = self['RawEventsTrialID'].astype('int32') - 1

        self['nTrials'] = int(self['RawEventsTrialID'][-1]) + 1  # Since index starts at zero

        # -- Fix the length of some arrays (FIX THIS, SHOULD NOT BE HARDCODED) --
        self['WithdrawalFromProbeOnset'] = self['WithdrawalFromProbeOnset'][:self['nTrials']]

    def mask_first_trial(self):
        eventsFirstTrial = np.flatnonzero(self['RawEventsTrialID']==1)
        self['RawEvents'][eventsFirstTrial,:] = 0

    def populate_dict(self):
        '''Load all data as dictionary key:value pairs'''
        for (key,value) in self.datanode._v_children.iteritems():
            if value.shape[0]>1:
                self[key] = value.read().transpose()
            else:
                self[key] = value.read()[0]

    def check_version(self):
        '''Verify that the file has the same version as this module'''
        if self.formatversion != FORMAT_VERSION_HDF5_MATLAB:
            errstr = 'Data format version (%s) is not %s.'%(self.formatversion,
                                                            FORMAT_VERSION_HDF5_MATLAB)
            raise TypeError(errstr)

    def find_stateID_each_trial(self,stateName):
        '''Returns the ID of a state on each trial. VERY INEFFICIENT.
        It returns -1 if the state did not exist for a given trial'''
        nTrials = self['nTrials']
        stringLength = self['StateLabelName'].shape[1]
        paddedName = string.ljust(stateName,stringLength)
        stateNameAsInt = np.array(struct.unpack('B'*stringLength,paddedName),dtype='uint8')
        stateID = np.empty(nTrials,dtype='int32')
        for indtrial in range(nTrials):
            limlo = np.searchsorted(self['StateLabelTrialID'],indtrial,side='left')
            limhi = np.searchsorted(self['StateLabelTrialID'],indtrial,side='right')
            thisTrialSlice = slice(limlo,limhi)
            thisTrialStateNames = self['StateLabelName'][thisTrialSlice,:]
            diffWithstateName = np.sum(thisTrialStateNames-stateNameAsInt, axis=1)
            indFirstStateName = np.flatnonzero(diffWithstateName==0)
            if indFirstStateName:
                thisTrialStateValues = self['StateLabelValue'][thisTrialSlice]
                stateID[indtrial] = thisTrialStateValues[indFirstStateName]
            else:
                stateID[indtrial] = -1
        return stateID

    def find_stateID(self,stateName):
        '''Returns the ID of a state (from the last trial).'''
        stringLength = self['StateLabelName'].shape[1]
        paddedName = string.ljust(stateName,stringLength)
        stateNameAsInt = np.array(struct.unpack('B'*stringLength,paddedName),dtype='uint8')
        indtrial = 1  # Not the empty one, but the first valid
        limlo = np.searchsorted(self['StateLabelTrialID'],indtrial,side='left')
        limhi = np.searchsorted(self['StateLabelTrialID'],indtrial,side='right')
        thisTrialSlice = slice(limlo,limhi)
        thisTrialStateNames = self['StateLabelName'][thisTrialSlice,:]
        diffWithstateName = np.sum(thisTrialStateNames-stateNameAsInt, axis=1)
        indFirstStateName = np.flatnonzero(diffWithstateName==0)
        if len(indFirstStateName):
            thisTrialStateValues = self['StateLabelValue'][thisTrialSlice]
            return thisTrialStateValues[indFirstStateName]
        else:
            raise ValueError('State %s does not exist.'%stateName)

    def time_start_each_trial(self,firstStateID):
        '''Returns the time of the transition into firstState for each trial.
        It assumes that the first state has the same ID on all trials.
        '''
        stateZeroID = 0 ##### WARNING!!! HARDCODED #####
        firstStateEventInds = np.logical_and(self['RawEvents'][:,0]==stateZeroID,
                                                self['RawEvents'][:,3]==firstStateID)
        firstStateEventTimes = self['RawEvents'][firstStateEventInds,2]
        return firstStateEventTimes

    def time_of_state_transition(self,prevStateID,nextStateID):
        '''Returns the time of the transition from prevStateID to nextStateID on each trial.
        It assumes that the states have the same ID on all trials.
        If the transition did not occur, it returns NaN for that trial.
        '''
        transitionEventTimes = np.empty(self['nTrials'])
        transitionEventTimes.fill(np.nan)
        transitionEventInds = np.logical_and(self['RawEvents'][:,0]==prevStateID,
                                             self['RawEvents'][:,3]==nextStateID)
        transitionTrialInd = self['RawEventsTrialID'][transitionEventInds]
        transitionEventTimes[transitionTrialInd] = self['RawEvents'][transitionEventInds,2]
        return transitionEventTimes

    def time_of_event(self,stateID,actionID):
        '''Returns the time of a given action within a given state (on each trial)
        It ignores all actions after the first one detected.
        It assumes that the states and actions have the same ID on all trials.
        If the event did not occur, it returns NaN for that trial.
        '''
        eventTimes = np.empty(self['nTrials'])
        eventTimes.fill(np.nan)
        eventInds = np.logical_and(self['RawEvents'][:,0]==stateID,
                                   self['RawEvents'][:,1]==actionID)
        transitionTrialInd = self['RawEventsTrialID'][eventInds]
        # Hack to ensure the first action is used instead of the last one.
        # A better solution is to use 'unique' with 'return_index' but this is
        # not available in Numpy_1.3
        eventInds = np.flatnonzero(eventInds)[::-1]
        transitionTrialInd = transitionTrialInd[::-1]
        eventTimes[transitionTrialInd] = self['RawEvents'][eventInds,2]
        return eventTimes
    


'''
behavData['RawEvents'][:20,:]
behavData['StateLabelName'][:10]
behavData['StateLabelTrialID']==1

eventInds = [2,3,5,20290]
transitionTrialInd = behavData['RawEventsTrialID'][eventInds]
unique(transitionTrialInd,return_index=True)

rawnode = h5file.getNode('/SessionParams/RawEvents')
rawdata = rawnode.read()
rawdata.shape

/var/data/test/data_saja_reversal_santiago_saja064_20110204a.h5 (File) ''
Last modif.: 'Fri Feb  4 15:41:26 2011'
Object Tree: 
/ (RootGroup) ''
/DataFormatVersion (Array(1,)) ''
/SessionParams (Group) ''
/SessionParams/ActionLabels (Table(1,)) ''
/SessionParams/AntiBiasMethod (CArray(1, 721), zlib(9)) ''
/SessionParams/AntiBiasMethodLabels (Table(1,)) ''
'''
