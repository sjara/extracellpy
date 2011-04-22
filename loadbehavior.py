#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Read behavior data'''

import numpy as np

__author__ = 'Santiago Jaramillo'
__version__ = '0.2'

import tables
import os
import string
import struct

def pp(rawev):
    '''Pretty print a matrix of raw events.
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

    def __init__(self,behavFile):
        h5file = tables.openFile(behavFile, mode = "r")
        self.fileName = behavFile
        self.formatversion = h5file.getNode('/DataFormatVersion').read()[0]
        self.check_version()
        self._datanode = h5file.getNode('/SessionParams') # Node is gone when closing file
        self.populate_dict()
        h5file.close()

        # -- Change first index to zero (the Python/C way) --
        self['StateLabelTrialID'] = self['StateLabelTrialID'].astype('int32') - 1
        self['RawEventsTrialID'] = self['RawEventsTrialID'].astype('int32') - 1

        self['nTrials'] = int(self['RawEventsTrialID'][-1]) + 1  # Since index starts at zero

    def mask_first_trial(self):
        eventsFirstTrial = np.flatnonzero(self['RawEventsTrialID']==1)
        self['RawEvents'][eventsFirstTrial,:] = 0

    def populate_dict(self):
        '''Load all data as dictionary paramName:value pairs'''
        for (paramName,paramValue) in self._datanode._v_children.iteritems():
            #print type(paramValue), paramName, paramValue.shape
            if isinstance(paramValue,tables.carray.CArray):
                if paramValue.shape[0]>1:
                    self[paramName] = paramValue.read().transpose()
                else:
                    self[paramName] = paramValue.read()[0]
            elif isinstance(paramValue,tables.array.Array):
                self[paramName] = paramValue.read()[0]
            elif isinstance(paramValue,tables.table.Table):
                dkeys   = paramValue.colnames
                dvalues = map(int,paramValue.read()[0])
                self[paramName] = dict(zip(dkeys,dvalues))
            else:
                raise TypeError('Unknown type: %s',paramName)
        #1/0

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
            #raise ValueError('State %s does not exist.'%stateName)
            print 'State %s does not exist.'%stateName
            return -1

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
        # FIXME: I'm changing the value of a parameter passed by reference
        # -- If input is a state name, find its ID --
        if isinstance(prevStateID,str):
            prevStateID = self.find_stateID(prevStateID)
        if isinstance(nextStateID,str):
            nextStateID = self.find_stateID(nextStateID)
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
        # FIXME: I'm changing the value of a parameter passed by reference
        # -- If input is a name, find its ID --
        if isinstance(stateID,str):
            stateID = self.find_stateID(stateID)
        if isinstance(actionID,str):
            actionID = self.find_actionID(actionID)
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

    def find_actionID(self,actionName):
        raise TypeError('This method is not implemented yet')

    def align_to_ephys(self,trialStartEphys):
        '''Find time of start of each trial according to the electrophysiology clock.'''
        # WARNING: this applies only to BControl with empty first trial 
        trialStartTime = self.trialStartTime[1:]
        # -- Find InterTrialIntervals --
        behavITI=np.diff(trialStartTime)
        ephysITI=np.diff(trialStartEphys)
        behavITI=np.hstack((0,behavITI))
        ephysITI=np.hstack((0,ephysITI))
        # -- Find closes ITI from behavior --
        bestInd = np.empty(behavITI.shape,dtype=np.intp)
        for indt,iti in enumerate(behavITI):
            bestInd[indt] = np.argmin(np.abs(iti-ephysITI))
        # -- Fix indexes for which ISI calculation yielded something off --
        for ind in range(1,bestInd.size-1):
            if abs(bestInd[ind]-bestInd[ind-1])>2:
                if bestInd[ind+1]-bestInd[ind-1]==2:
                    bestInd[ind] = bestInd[ind-1]+1
        #ephysTimeStamps = np.empty(trialStartTime.shape)
        #ephysTimeStamps.fill(np.NaN)
        #ephysTimeStamps= trialStartTimeNL[bestInd]
        self.trialStartTimeEphys = trialStartEphys[bestInd]
        # WARNING: this applies only to BControl with empty first trial 
        self.trialStartTimeEphys = np.hstack((self.trialStartTimeEphys[0],
                                              self.trialStartTimeEphys))

    def OTHER_align_to_ephys(self,trialStartEphys,plot=0):
        '''Find time of start of each trial according to the electrophysiology clock.'''
        self.trialStartTimeEphys = trialStartEphys[:-1]
        # WARNING: this applies only to BControl with empty first trial 
        self.trialStartTimeEphys = np.hstack((self.trialStartTimeEphys[0],
                                              self.trialStartTimeEphys))

    def check_clock_drift(self):
        '''Plot comparison between behavior and electrophysiology clocks.'''
        import pylab as p
        behavTS = self.trialStartTime
        ephysTS = self.trialStartTimeEphys
        p.clf()
        p.subplot(2,1,1)
        p.hold(True)
        pB = p.plot(range(behavTS.size),behavTS-behavTS[1],'o',mec='b',mfc='none')
        pE = p.plot(range(ephysTS.size),ephysTS-ephysTS[1],'.r')
        p.hold(False)
        p.ylabel('Time w.r.t. first event (sec)')
        p.legend((pB,pE),('Behavior time','Neuralynx time'),numpoints=1,loc='upper left')
        p.subplot(2,1,2)
        drift = (ephysTS-ephysTS[1])-(behavTS-behavTS[1])
        p.plot(range(behavTS.size),drift,'o',mec='b',mfc='none')
        p.ylabel('Drift between the two clocks (sec)')
        p.xlabel('Trial')
        p.show()


class ReversalBehaviorData(BehaviorData):
    '''This class inherits BehaviorData and adds methods specific to reversal protocol.'''
    def __init__(self,behavFileName):
        BehaviorData.__init__(self,behavFileName)
        self['HitHistoryLabels'] = {'Correct':1, 'Error':0, 'EarlyWithdrawal':-1,
                                     'ErrorNextCorr':-2, 'CorrectNextCorr':2, 'Direct':2,
                                     'TimeOut':-3}
        # -- Remove irrelevant values at the end of some arrays --
        arraysToFix = ['HitHistory','WithdrawalFromProbeOnset',
                       'RewardSideList','PreStimTime','CurrentBlock','TargetFreq']
        for arrayName in arraysToFix:
            self[arrayName] = self[arrayName][:self['nTrials']]
    def extract_event_times(self):
        # NOTE: ActionLabels does not have correct labels (Cin=1, Cout=2, ...)
        #       Maybe BControl's get_col_labels(current_assembler) is wrong.
        # -- Santiago's rig --
        centerPokeOutID = 1 ### WARNING!!! HARDCODED
        leftPokeInID    = 2 ### WARNING!!! HARDCODED
        rightPokeInID   = 4 ### WARNING!!! HARDCODED
        '''
        # -- Peter's rig --
        centerPokeOutID = 2 ### WARNING!!! HARDCODED
        leftPokeInID    = 16 ### WARNING!!! HARDCODED
        rightPokeInID   = 4 ### WARNING!!! HARDCODED
        '''
        self.trialStartTime = self.time_of_state_transition('state_0','send_trial_info')
        #self.trialStartTime = self.time_of_state_transition(16,'send_trial_info')
        self.targetOnsetTime = self.time_of_state_transition('delay_period','play_target')
        self.centerOutTime = self.time_of_event('wait_for_apoke',centerPokeOutID)
        self.leftInTime = self.time_of_event('wait_for_apoke',leftPokeInID)
        self.rightInTime = self.time_of_event('wait_for_apoke',rightPokeInID)
    def find_trials_each_type(self):
        self.highFreqs = self['CurrentBlock']==self['CurrentBlockLabels']['HighFreqs']
        self.lowFreqs = self['CurrentBlock']==self['CurrentBlockLabels']['LowFreqs']
        self.leftReward = self['RewardSideList']==self['RewardSideListLabels']['left']
        self.rightReward = self['RewardSideList']==self['RewardSideListLabels']['right']
        self.correct = self['HitHistory']==self['HitHistoryLabels']['Correct']
        self.error = self['HitHistory']==self['HitHistoryLabels']['Error']
        self.early = self['HitHistory']==self['HitHistoryLabels']['EarlyWithdrawal']
        #self. = self['']==self['Labels']['']
        self.leftChoice = np.logical_not(np.isnan(self.leftInTime))
        self.rightChoice = np.logical_not(np.isnan(self.rightInTime))
        self.sideInTime = np.copy(self.leftInTime)
        self.sideInTime[self.rightChoice] = self.rightInTime[self.rightChoice]
    def find_boundaries_each_block(self):
        # FIXME: for efficiency, check if lastTrialEachBlock already exists
        blockBoundaries = np.flatnonzero(np.diff(self['CurrentBlock']))
        self.lastTrialEachBlock = np.hstack((blockBoundaries,self['nTrials']))
        self.firstTrialEachBlock = np.hstack((0,self.lastTrialEachBlock[:-1]+1))
        self['nBlocks'] = len(self.lastTrialEachBlock)
    def find_trials_each_block(self):
        self.find_boundaries_each_block()
        self.trialsEachBlock = np.zeros((self['nTrials'],self['nBlocks']),dtype='bool')
        for block in range(self['nBlocks']):
            bSlice = slice(self.firstTrialEachBlock[block],self.lastTrialEachBlock[block]+1)
            self.trialsEachBlock[bSlice,block]=True
        

class TuningBehaviorData(BehaviorData):
    '''This class inherits BehaviorData and adds methods specific to tuningcurve protocol.'''
    def __init__(self,behavFileName):
        BehaviorData.__init__(self,behavFileName)
        arraysToFix = ['SoundFreq']
        for arrayName in arraysToFix:
            self[arrayName] = self[arrayName][:self['nTrials']]
    def extract_event_times(self):
        # NOTE: ActionLabels does not have correct labels (Cin=1, Cout=2, ...)
        #       Maybe BControl's get_col_labels(current_assembler) is wrong.
        centerPokeOutID = 1 ### WARNING!!! HARDCODED
        leftPokeInID    = 2 ### WARNING!!! HARDCODED
        rightPokeInID   = 4 ### WARNING!!! HARDCODED
        self.trialStartTime = self.time_of_state_transition('state_0','send_trial_info')
        self.targetOnsetTime = self.time_of_state_transition('continue_trial','play_sound')



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
