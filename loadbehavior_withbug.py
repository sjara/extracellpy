#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Read behavior data

USE THIS VERSION FOR DATA PRE 2011-03-19 
For these data files:
 'HitHistory','RewardSideList','WithdrawalFromProbeOnset':saved from 1st element
 'TargetFreq','CurrentBlock': saved from 2nd element
'''

import numpy as np

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

import tables
import os
import string
import struct
import loadbehavior

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


class ReversalBehaviorData(loadbehavior.ReversalBehaviorData):
    '''This class inherits BehaviorData and adds methods specific to reversal protocol.'''
    def __init__(self,behavFileName):
        loadbehavior.ReversalBehaviorData.__init__(self,behavFileName)
        # -- Fixes for buggy data --
        self.remove_trial_from_raw(1)
        self.fix_bug()
    def remove_trial_from_raw(self,trialToRemove):
        '''Careful, this does not remove the data on all other vectors.'''
        eventsToKeep = np.logical_not(self['RawEventsTrialID']==trialToRemove)
        self['RawEvents'] = self['RawEvents'][eventsToKeep,:]
        trialIDtoChange = self['RawEventsTrialID']>trialToRemove
        self['RawEventsTrialID'][trialIDtoChange] = self['RawEventsTrialID'][trialIDtoChange]-1
        self['RawEventsTrialID'] = self['RawEventsTrialID'][eventsToKeep,:]
        self['nTrials'] = self['nTrials']-1
        # FIXME: Do I need to change 'StateLabelTrialID'?
        #behavData['StateLabelTrialID']
    def fix_bug(self):
        '''This method fixes the misalignment between vectors on data pre 2011-03-19'''
        # Remove first trial (because other variables do not have it)
        self['WithdrawalFromProbeOnset'] = self['WithdrawalFromProbeOnset'][1:]
        self['HitHistory'] = self['HitHistory'][1:]
        self['RewardSideList'] = self['RewardSideList'][1:]
        # Remove last trial (because the trials was not completed)
        self['TargetFreq'] = self['TargetFreq'][0:-1]
        self['CurrentBlock'] = self['CurrentBlock'][0:-1]
        self['ProbeDuration'] = self['ProbeDuration'][0:-1]
        self['DelayToTarget'] = self['DelayToTarget'][0:-1]
        self['TargetDuration'] = self['TargetDuration'][0:-1]
    def extract_event_times(self):
        # NOTE: ActionLabels does not have correct labels (Cin=1, Cout=2, ...)
        #       Maybe BControl's get_col_labels(current_assembler) is wrong.
        #       so instead, here are the definitions for each rig.
        rig = self['HostName']
        if rig in ['lermontov','peter']:
            # -- Peter's rig --
            centerPokeInID  = 1 ### WARNING!!! HARDCODED
            centerPokeOutID = 2 ### WARNING!!! HARDCODED
            leftPokeInID    = 4 ### WARNING!!! HARDCODED (What I thought before 2011-09-25)
            rightPokeInID   = 16 ### WARNING!!! HARDCODED (What I thought before 2011-09-25)
            #leftPokeInID    = 16 ### WARNING!!! HARDCODED (What I thought before 2011-09-25)
            #rightPokeInID   = 4 ### WARNING!!! HARDCODED (What I thought before 2011-09-25)
        elif rig in ['cnmcx','cnmc5','cnmc6','cnmc7','cnmc8','cnmc10','cnmc11','cnmc12']:
            # -- Most rigs (and old Linux state machine) --
            centerPokeInID = 0 ### WARNING!!! HARDCODED
            centerPokeOutID = 1 ### WARNING!!! HARDCODED
            leftPokeInID    = 4 ### WARNING!!! HARDCODED
            rightPokeInID   = 16 ### WARNING!!! HARDCODED
        elif rig in ['cnmc9']:
            # -- Santiago's rig (with newer Linux state machine) --
            # Except that back in the day, it was not different that the others
            centerPokeInID = 0 ### WARNING!!! HARDCODED
            centerPokeOutID = 1 ### WARNING!!! HARDCODED
            leftPokeInID    = 4#2 ### WARNING!!! HARDCODED
            rightPokeInID   = 16# ### WARNING!!! HARDCODED
        else:
            raise ValueError('Rig number not defined')
        self.trialStartTime = self.time_of_state_transition('state_0','send_trial_info')
        #self.trialStartTime = self.time_of_state_transition(16,'send_trial_info')
        self.targetOnsetTime = self.time_of_state_transition('delay_period','play_target')
        self.centerInTime = self.time_of_event('wait_for_cpoke',centerPokeInID)
        self.centerOutTime = self.time_of_event('wait_for_apoke',centerPokeOutID)
        self.leftInTime = self.time_of_event('wait_for_apoke',leftPokeInID)
        self.rightInTime = self.time_of_event('wait_for_apoke',rightPokeInID)



class OLD_ReversalBehaviorData(loadbehavior.BehaviorData):
    '''This class inherits BehaviorData and adds methods specific to reversal protocol.'''
    def __init__(self,behavFileName):
        BehaviorData.__init__(self,behavFileName)
        self['HitHistoryLabels'] = {'Correct':1, 'Error':0, 'EarlyWithdrawal':-1,
                                     'ErrorNextCorr':-2, 'CorrectNextCorr':2, 'Direct':2,
                                     'TimeOut':-3}
        # -- Remove irrelevant values at the end of some arrays --
        arraysToFix = ['HitHistory','WithdrawalFromProbeOnset',
                       'RewardSideList','PreStimTime']
        for arrayName in arraysToFix:
            self[arrayName] = self[arrayName][:self['nTrials']]
        # -- Fixes for buggy data --
        self.remove_trial_from_raw(1)
        self.fix_bug()
    def fix_bug(self):
        '''This method fixes the misalignment between vectors on data pre 2011-03-19'''
        # Remove first trial (because other variables do not have it)
        self['WithdrawalFromProbeOnset'] = self['WithdrawalFromProbeOnset'][1:]
        self['HitHistory'] = self['HitHistory'][1:]
        self['RewardSideList'] = self['RewardSideList'][1:]
        # Remove last trial (because the trials was not completed)
        self['TargetFreq'] = self['TargetFreq'][0:-1]
        self['CurrentBlock'] = self['CurrentBlock'][0:-1]
        self['ProbeDuration'] = self['ProbeDuration'][0:-1]
        self['DelayToTarget'] = self['DelayToTarget'][0:-1]
        self['TargetDuration'] = self['TargetDuration'][0:-1]
    def extract_event_times(self):
        '''
        # NOTE: ActionLabels does not have correct labels (Cin=1, Cout=2, ...)
        #       Maybe BControl's get_col_labels(current_assembler) is wrong.
        centerPokeOutID = 1 ### WARNING!!! HARDCODED
        leftPokeInID    = 2 ### WARNING!!! HARDCODED
        rightPokeInID   = 4 ### WARNING!!! HARDCODED
        '''
        # NOTE: ActionLabels does not have correct labels (Cin=1, Cout=2, ...)
        #       Maybe BControl's get_col_labels(current_assembler) is wrong.
        #       so instead, here are the definitions for each rig.
        rig = self['HostName']
        if rig in ['lermontov','peter']:
            # -- Peter's rig --
            centerPokeInID  = 1 ### WARNING!!! HARDCODED
            centerPokeOutID = 2 ### WARNING!!! HARDCODED
            leftPokeInID    = 4 ### WARNING!!! HARDCODED (What I thought before 2011-09-25)
            rightPokeInID   = 16 ### WARNING!!! HARDCODED (What I thought before 2011-09-25)
            #leftPokeInID    = 16 ### WARNING!!! HARDCODED (What I thought before 2011-09-25)
            #rightPokeInID   = 4 ### WARNING!!! HARDCODED (What I thought before 2011-09-25)
        elif rig in ['cnmcx','cnmc5','cnmc6','cnmc7','cnmc8','cnmc10','cnmc11','cnmc12']:
            # -- Most rigs (and old Linux state machine) --
            centerPokeInID = 0 ### WARNING!!! HARDCODED
            centerPokeOutID = 1 ### WARNING!!! HARDCODED
            leftPokeInID    = 4 ### WARNING!!! HARDCODED
            rightPokeInID   = 16 ### WARNING!!! HARDCODED
        elif rig in ['cnmc9']:
            # -- Santiago's rig (with newer Linux state machine) --
            centerPokeInID = 0 ### WARNING!!! HARDCODED
            centerPokeOutID = 1 ### WARNING!!! HARDCODED
            leftPokeInID    = 4#2 ### WARNING!!! HARDCODED
            rightPokeInID   = 16#4 ### WARNING!!! HARDCODED
        else:
            raise ValueError('Rig number not defined')
        self.trialStartTime = self.time_of_state_transition('state_0','send_trial_info')
        #self.trialStartTime = self.time_of_state_transition(16,'send_trial_info')
        self.targetOnsetTime = self.time_of_state_transition('delay_period','play_target')
        self.centerInTime = self.time_of_event('wait_for_cpoke',centerPokeInID)
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
