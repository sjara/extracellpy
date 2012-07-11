#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Functions for analysis of electrophysiology/behavior session.
'''

from extracellpy import settings
reload(settings) # Force reload
from extracellpy import loadneuralynx
reload(loadneuralynx) # Force reload
from extracellpy import spikesanalysis
from extracellpy import spikesdetection
from extracellpy import loadbehavior
from extracellpy import colorpalette as cp
import os, sys
import numpy as np

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

BEHAVIORPATH = settings.BEHAVIOR_PATH
EPHYSPATH = settings.EPHYS_PATH
EXTRACTED_SPIKES_PATH = settings.EXTRACTED_SPIKES_PATH

bitTRIALIND_DEFAULT = 3      # TrialIndicator (bitID starting from 0)
bitPHOTOSTIMIND = 4  # PhotoStimIndicator (bitID starting from 0)
bitTARGETIND = 5     # TargetIndicator (bitID starting from 0)

npAND = np.logical_and
npOR = np.logical_or
npNOT = np.logical_not


def load_behavior_session(animalName,behavSession):
    behavDataDir = os.path.join(BEHAVIORPATH,'%s/'%(animalName))
    behavFileName = 'data_saja_reversal_santiago_%s_%s.h5'%(animalName,behavSession)
    behavFile = os.path.join(behavDataDir,behavFileName)
    behavData = loadbehavior.ReversalBehaviorData(behavFile)
    behavData.extract_event_times()
    behavData.find_trials_each_type()
    return behavData

def get_behavior_filename(animalName,behavSession):
    behavDataDir = os.path.join(BEHAVIORPATH,'%s/'%(animalName))
    behavFileName = 'data_saja_reversal_santiago_%s_%s.h5'%(animalName,behavSession)
    behavFile = os.path.join(behavDataDir,behavFileName)
    return behavFile

def load_behavior_one_cell(oneCell):
    animalName = oneCell.animalName
    ephysSession = oneCell.ephysSession
    behavSession = oneCell.behavSession
    #behavDataDir = '/var/data/BControlData/Data/santiago/%s/'%(animalName)
    behavDataDir = os.path.join(BEHAVIORPATH,'%s/'%(animalName))
    behavFileName = 'data_saja_reversal_santiago_%s_%s.h5'%(animalName,behavSession)
    behavFile = os.path.join(behavDataDir,behavFileName)
    behavData = loadbehavior.ReversalBehaviorData(behavFile)
    behavData.extract_event_times()
    behavData.find_trials_each_type()
    #behavData.align_to_ephys(trialStartTimeNL)
    return behavData

def load_lfp_reversal(oneLFP,prevLFP=None,prevData=None,bitTRIALIND=bitTRIALIND_DEFAULT):
    '''Load behavior and LFP data.
    FIXME: this function replicates many things from load_cell_reversal, they should be unified
    '''
    animalName = oneLFP.animalName
    ephysSession = oneLFP.ephysSession
    behavSession = oneLFP.behavSession
    electrode = oneLFP.electrode
    # -- Load events from Neuralynx --
    if not prevLFP or oneCell.ephysSession!=prevLFP.ephysSession:
        dataDir = os.path.join(EPHYSPATH,'%s/%s/'%(animalName,ephysSession))
        eventsFile = os.path.join(dataDir,'Events.nev')
        events = loadneuralynx.DataEvents(eventsFile)
        trialEvents = (events.valueTTL & (1<<bitTRIALIND)) != 0
        trialStartTimeNL = 1e-6*events.timestamps[trialEvents]
        #timeOfFirstTrial = trialStartTimeNL[0]  # First non-empty trial
        #print '******** loading ephys **********'
        ### NOTE: trialStartTimeNL does not contain the first empty behavior trial.
        ###       First behavior trial is removed in behavData.align_to_ephys()
    else:
        trialEvents = prevData['trialEvents']
    # -- Load events from behavior --
    if not prevLFP or oneCell.behavSession!=prevLFP.behavSession:
        behavDataDir = os.path.join(BEHAVIORPATH,'%s/'%(animalName))
        behavFileName = 'data_saja_reversal_santiago_%s_%s.h5'%(animalName,behavSession)
        behavFile = os.path.join(behavDataDir,behavFileName)
        behavData = loadbehavior.ReversalBehaviorData(behavFile)
        behavData.extract_event_times()
        behavData.find_trials_each_type()
        behavData.align_to_ephys(trialStartTimeNL)
        #print '******** loading behavior **********'
    else:
        behavData = prevData['behavData']
    # -- Load LFP --
    if not prevLFP or oneCell.tetrode!=prevLFP.tetrode:
        contFile = os.path.join(dataDir,'CSC%d.ncs'%electrode)
        dataLFP = loadneuralynx.DataCont(contFile)
    else:
        dataLFP = prevData['dataLFP']
    return (behavData,trialEvents,dataLFP)
     

def load_cell_reversal(oneCell,prevCell=None,prevData=None,
                       bitTRIALIND=bitTRIALIND_DEFAULT,clusterFileSuffix='1'):
    '''Load behavior and spikes data.
       returns (behavData,trialEvents,dataTT,spikeInds)
    '''
    #ephysSession = '2011-05-10_17-05-22'
    #behavSession = '20110510a'
    animalName = oneCell.animalName
    ephysSession = oneCell.ephysSession
    behavSession = oneCell.behavSession
    tetrode = oneCell.tetrode
    cluster = oneCell.cluster

    # -- Load events from Neuralynx --
    if not prevCell or oneCell.ephysSession!=prevCell.ephysSession:
        dataDir = os.path.join(EPHYSPATH,'%s/%s/'%(animalName,ephysSession))
        clustersDir = os.path.join(EPHYSPATH,'%s/%s_kk/'%(animalName,ephysSession))
        eventsFile = os.path.join(dataDir,'Events.nev')
        events = loadneuralynx.DataEvents(eventsFile)
        trialEvents = (events.valueTTL & (1<<bitTRIALIND)) != 0
        trialStartTimeNL = 1e-6*events.timestamps[trialEvents]
        #timeOfFirstTrial = trialStartTimeNL[0]  # First non-empty trial
        #print '******** loading ephys **********'
        ### NOTE: trialStartTimeNL does not contain the first empty behavior trial.
        ###       First behavior trial is removed in behavData.align_to_ephys()
    else:
        trialEvents = prevData['trialEvents']
    # -- Load events from behavior --
    if not prevCell or oneCell.behavSession!=prevCell.behavSession:
        behavDataDir = os.path.join(BEHAVIORPATH,'%s/'%(animalName))
        behavFileName = 'data_saja_reversal_santiago_%s_%s.h5'%(animalName,behavSession)
        behavFile = os.path.join(behavDataDir,behavFileName)
        behavData = loadbehavior.ReversalBehaviorData(behavFile)
        behavData.extract_event_times()
        behavData.find_trials_each_type()
        behavData.align_to_ephys(trialStartTimeNL)
        #print '******** loading behavior **********'
    else:
        behavData = prevData['behavData']
    # -- Load spikes --
    if not prevCell or oneCell.tetrode!=prevCell.tetrode:
        tetrodeFile = os.path.join(dataDir,'TT%d.ntt'%tetrode)
        dataTT = loadneuralynx.DataTetrode(tetrodeFile)
        dataTT.timestamps = dataTT.timestamps.astype(np.float64)*1e-6  # in sec
    else:
        dataTT = prevData['dataTT']
    # -- Load clusters --
    clustersFile = os.path.join(clustersDir,'TT%d.clu.%s'%(tetrode,clusterFileSuffix))
    dataTT.set_clusters(clustersFile)
    spikeInds = np.flatnonzero(dataTT.clusters==cluster)

    return (behavData,trialEvents,dataTT,spikeInds)


def load_cell_reversal_fromcont(oneCell,prevCell=None,prevData=None,bitTRIALIND=bitTRIALIND_DEFAULT):
    '''Load behavior and spikes data. Spikes data in this case comes from detected spikes
    from continuous data.
    '''
    #ephysSession = '2011-05-10_17-05-22'
    #behavSession = '20110510a'
    animalName = oneCell.animalName
    ephysSession = oneCell.ephysSession
    behavSession = oneCell.behavSession
    tetrode = oneCell.tetrode
    cluster = oneCell.cluster

    # -- Load events from Neuralynx --
    if not prevCell or oneCell.ephysSession!=prevCell.ephysSession:
        dataDir = os.path.join(EPHYSPATH,'%s/%s/'%(animalName,ephysSession))
        #clustersDir = os.path.join(EPHYSPATH,'%s/%s_kk/'%(animalName,ephysSession))
        eventsFile = os.path.join(dataDir,'Events.nev')
        events = loadneuralynx.DataEvents(eventsFile)
        trialEvents = (events.valueTTL & (1<<bitTRIALIND)) != 0
        trialStartTimeNL = 1e-6*events.timestamps[trialEvents]
        #print '******** loading ephys **********'
        ### NOTE: trialStartTimeNL does not contain the first empty behavior trial.
        ###       First behavior trial is removed in behavData.align_to_ephys()
    else:
        trialEvents = prevData['trialEvents']
    # -- Load events from behavior --
    if not prevCell or oneCell.behavSession!=prevCell.behavSession:
        behavDataDir = os.path.join(BEHAVIORPATH,'%s/'%(animalName))
        behavFileName = 'data_saja_reversal_santiago_%s_%s.h5'%(animalName,behavSession)
        behavFile = os.path.join(behavDataDir,behavFileName)
        behavData = loadbehavior.ReversalBehaviorData(behavFile)
        behavData.extract_event_times()
        behavData.find_trials_each_type()
        behavData.align_to_ephys(trialStartTimeNL)
        #print '******** loading behavior **********'
    # -- Load spikes --
    if not prevCell or oneCell.tetrode!=prevCell.tetrode:
        spikesDataDir = EXTRACTED_SPIKES_PATH%animalName
        if cluster is None:
            spikesFileName = '%s_%s_e%02d_spikes.h5'%(animalName,ephysSession,tetrode)
        else:
            spikesFileName = '%s_%s_e%02d_c%d_spikes.h5'%(animalName,ephysSession,tetrode,cluster)
        spikesFileFull = os.path.join(spikesDataDir,spikesFileName)
        dataTT = spikesdetection.ExtractedSpikes(spikesFileFull)
        dataTT.timestamps = dataTT.timestamps.astype(np.float64)*1e-6  # in sec
    else:
        dataTT = prevData['dataTT']
    '''
    # -- Load clusters --
    clustersFile = os.path.join(clustersDir,'TT%d.clu.1'%tetrode)
    dataTT.set_clusters(clustersFile)
    spikeInds = np.flatnonzero(dataTT.clusters==cluster)
    '''
    spikeInds = None
    return (behavData,trialEvents,dataTT,spikeInds)

def load_mu_reversal(oneSite,prevSite=None,prevData=None,bitTRIALIND=bitTRIALIND_DEFAULT):
    '''Load behavior and spikes data.
    '''
    import extrafuncs
    animalName = oneSite.animalName
    ephysSession = oneSite.ephysSession
    behavSession = oneSite.behavSession
    tetrode = oneSite.tetrode
    clusters = oneSite.clusters

    # -- Load events from Neuralynx --
    if not prevSite or oneSite.ephysSession!=prevSite.ephysSession:
        dataDir = os.path.join(EPHYSPATH,'%s/%s/'%(animalName,ephysSession))
        clustersDir = os.path.join(EPHYSPATH,'%s/%s_kk/'%(animalName,ephysSession))
        eventsFile = os.path.join(dataDir,'Events.nev')
        events = loadneuralynx.DataEvents(eventsFile)
        trialEvents = (events.valueTTL & (1<<bitTRIALIND)) != 0
        trialStartTimeNL = 1e-6*events.timestamps[trialEvents]
        #timeOfFirstTrial = trialStartTimeNL[0]  # First non-empty trial
        #print '******** loading ephys **********'
        ### NOTE: trialStartTimeNL does not contain the first empty behavior trial.
        ###       First behavior trial is removed in behavData.align_to_ephys()
    else:
        trialEvents = prevData['trialEvents']
    # -- Load events from behavior --
    if not prevSite or oneSite.behavSession!=prevSite.behavSession:
        behavDataDir = os.path.join(BEHAVIORPATH,'%s/'%(animalName))
        behavFileName = 'data_saja_reversal_santiago_%s_%s.h5'%(animalName,behavSession)
        behavFile = os.path.join(behavDataDir,behavFileName)
        behavData = loadbehavior.ReversalBehaviorData(behavFile)
        behavData.extract_event_times()
        behavData.find_trials_each_type()
        behavData.align_to_ephys(trialStartTimeNL)
        #print '******** loading behavior **********'
    else:
        behavData = prevData['behavData']
    # -- Load spikes --
    if not prevSite or oneSite.tetrode!=prevSite.tetrode:
        tetrodeFile = os.path.join(dataDir,'TT%d.ntt'%tetrode)
        dataTT = loadneuralynx.DataTetrode(tetrodeFile)
        dataTT.timestamps = dataTT.timestamps.astype(np.float64)*1e-6  # in sec
    else:
        dataTT = prevData['dataTT']
    # -- Load clusters --
    if clusters:
        clustersFile = os.path.join(clustersDir,'TT%d.clu.1'%tetrode)
        dataTT.set_clusters(clustersFile)
        #spikeInds = np.flatnonzero(dataTT.clusters==cluster)
        spikeInds = np.flatnonzero(extrafuncs.ismember(dataTT.clusters,clusters))
    else:
        spikeInds = slice(len(dataTT.timestamps))
    '''
    spikeInds = []
    '''
    return (behavData,trialEvents,dataTT,spikeInds)


def load_cell_tuningcurve(oneCell,prevCell=None,prevData=None,bitTRIALIND=bitTRIALIND_DEFAULT):
    '''Load behavior and spikes data for cell recorded under saja_tuning protocol.
    '''
    animalName = oneCell.animalName
    ephysSession = oneCell.ephysSession
    behavSession = oneCell.behavSession
    tetrode = oneCell.tetrode
    cluster = oneCell.cluster

    # -- Load events from Neuralynx --
    if not prevCell or oneCell.ephysSession!=prevCell.ephysSession:
        dataDir = os.path.join(EPHYSPATH,'%s/%s/'%(animalName,ephysSession))
        clustersDir = os.path.join(EPHYSPATH,'%s/%s_kk/'%(animalName,ephysSession))
        eventsFile = os.path.join(dataDir,'Events.nev')
        events = loadneuralynx.DataEvents(eventsFile)
        eventsTime={}
        (trialOnset,trialOffset) = events.find_bit_changes(bitTRIALIND)
        (targetOnset,targetOffset) = events.find_bit_changes(bitTARGETIND)
        (photoOnset,photoOffset) = events.find_bit_changes(bitPHOTOSTIMIND)
        eventsTime['trial'] = 1e-6*events.timestamps[trialOnset]
        eventsTime['targetAlone'] = 1e-6*events.timestamps[npAND(targetOnset,npNOT(photoOnset))]
        eventsTime['targetAndPhotostim'] = 1e-6*events.timestamps[npAND(targetOnset,photoOnset)]
    else:
        trialEvents = prevData['trialEvents']
    '''
    # -- Load events from behavior --
    if not prevCell or oneCell.behavSession!=prevCell.behavSession:
        behavDataDir = '/var/data/BControlData/Data/santiago/%s/'%(animalName)
        behavFileName = 'data_saja_reversal_santiago_%s_%s.h5'%(animalName,behavSession)
        behavFile = os.path.join(behavDataDir,behavFileName)
        behavData = loadbehavior.ReversalBehaviorData(behavFile)
        behavData.extract_event_times()
        behavData.find_trials_each_type()
        behavData.align_to_ephys(trialStartTimeNL)
        #print '******** loading behavior **********'
    else:
        behavData = prevData['behavData']
    '''
    behavData = []
    # -- Load spikes --
    if not prevCell or oneCell.tetrode!=prevCell.tetrode:
        tetrodeFile = os.path.join(dataDir,'TT%d.ntt'%tetrode)
        dataTT = loadneuralynx.DataTetrode(tetrodeFile)
        dataTT.timestamps = dataTT.timestamps.astype(np.float64)*1e-6  # in sec
    else:
        dataTT = prevData['dataTT']
    # -- Load clusters --
    clustersFile = os.path.join(clustersDir,'TT%d.clu.1'%tetrode)
    dataTT.set_clusters(clustersFile)
    spikeInds = np.flatnonzero(dataTT.clusters==cluster)

    return (behavData,eventsTime,dataTT,spikeInds)


def trials_by_condition(behavData,CONDCASE=1,sortby=np.empty(0),outcome='all',selected=None):
    ''' 
    (trialsEachCond,condInfo) = trials_by_condition(behavData,1)
    '''
    if selected is None:
        selected = np.ones(len(behavData.correct),dtype=bool)
    if CONDCASE==1:
        trialsEachCondLabels = {'LowFreq:LeftReward':0, 'MidFreq:RightReward':1,
                                'MidFreq:LeftReward':2, 'HighFreq:RightReward':3}
        colorEachCond = [cp.TangoPalette['SkyBlue1'], cp.TangoPalette['Orange2'],
                         cp.TangoPalette['ScarletRed1'], cp.TangoPalette['Chameleon3']] #'Plum1'
        # FIXME: finish implementing this properly
        if outcome=='correct':
            g1 = np.flatnonzero(behavData.lowFreqs & behavData.leftReward & behavData.correct & selected)
            g2 = np.flatnonzero(behavData.lowFreqs & behavData.rightReward & behavData.correct & selected)
            g3 = np.flatnonzero(behavData.highFreqs & behavData.leftReward & behavData.correct & selected)
            g4 = np.flatnonzero(behavData.highFreqs & behavData.rightReward & behavData.correct & selected)
        elif outcome=='error':
            print 'WARNING! this needs fixing. It does not exclude bad trials'
            g1 = np.flatnonzero(npAND(npAND(behavData.lowFreqs,behavData.leftReward),behavData.error))
            g2 = np.flatnonzero(npAND(npAND(behavData.lowFreqs,behavData.rightReward),behavData.error))
            g3 = np.flatnonzero(npAND(npAND(behavData.highFreqs,behavData.leftReward),behavData.error))
            g4 = np.flatnonzero(npAND(npAND(behavData.highFreqs,behavData.rightReward),behavData.error))
        elif outcome=='early':
            print 'WARNING! this needs fixing. It does not exclude bad trials'
            g1 = np.flatnonzero(npAND(npAND(behavData.lowFreqs,behavData.leftReward),behavData.early))
            g2 = np.flatnonzero(npAND(npAND(behavData.lowFreqs,behavData.rightReward),behavData.early))
            g3 = np.flatnonzero(npAND(npAND(behavData.highFreqs,behavData.leftReward),behavData.early))
            g4 = np.flatnonzero(npAND(npAND(behavData.highFreqs,behavData.rightReward),behavData.early))
        elif outcome=='valid':
            print 'WARNING! this needs fixing. It does not exclude bad trials'
            g1 = np.flatnonzero(npAND(npAND(behavData.lowFreqs,behavData.leftReward),npNOT(behavData.early)))
            g2 = np.flatnonzero(npAND(npAND(behavData.lowFreqs,behavData.rightReward),npNOT(behavData.early)))
            g3 = np.flatnonzero(npAND(npAND(behavData.highFreqs,behavData.leftReward),npNOT(behavData.early)))
            g4 = np.flatnonzero(npAND(npAND(behavData.highFreqs,behavData.rightReward),npNOT(behavData.early)))
        elif outcome=='correctness':
            trialsEachCondLabels = {'MidFreq:RightReward':0, 'MidFreq:LeftReward':1,
                                    'MidFreq:RightRewardError':2, 'MidFreq:LeftRewardError':3}
            #colorEachCond = [cp.TangoPalette['Orange1'], cp.TangoPalette['ScarletRed1'],
            #                 cp.TangoPalette['Orange3'], cp.TangoPalette['ScarletRed3']]
            colorEachCond = [cp.TangoPalette['Orange1'], 'y',
                             cp.TangoPalette['ScarletRed1'],'m']
            g1 = np.flatnonzero(behavData.lowFreqs & behavData.rightReward & behavData.correct & selected)
            g2 = np.flatnonzero(behavData.lowFreqs & behavData.rightReward & behavData.error & selected)
            g3 = np.flatnonzero(behavData.highFreqs & behavData.leftReward & behavData.correct & selected)
            g4 = np.flatnonzero(behavData.highFreqs & behavData.leftReward & behavData.error & selected)
        elif outcome=='correctPerBlock':
            behavData.find_trials_each_block() # BAD CODING: changing object obscurely
            g1part = behavData.lowFreqs & behavData.leftReward & behavData.correct & selected
            g2part = behavData.lowFreqs & behavData.rightReward & behavData.correct & selected
            g3part = behavData.highFreqs & behavData.leftReward & behavData.correct & selected
            g4part = behavData.highFreqs & behavData.rightReward & behavData.correct & selected
            if np.mean(behavData.lowFreqs[:10])>np.mean(behavData.highFreqs[:10]):
                trialsEachCondLabels = ['14kHz (R)','14kHz (L)','14kHz (R)','14kHz (L)']
                colorEachCond = [cp.TangoPalette['Orange2'], cp.TangoPalette['ScarletRed1'],
                                 cp.TangoPalette['Orange2'], cp.TangoPalette['ScarletRed1']]
                #g1part = npAND(npAND(behavData.lowFreqs,behavData.leftReward),behavData.correct)
                #g2part = npAND(npAND(behavData.lowFreqs,behavData.rightReward),behavData.correct)
                #g3part = npAND(npAND(behavData.highFreqs,behavData.leftReward),behavData.correct)
                #g4part = npAND(npAND(behavData.highFreqs,behavData.rightReward),behavData.correct)
                g1 = np.flatnonzero(npAND(g2part,behavData.trialsEachBlock[:,1]))
                g3 = np.flatnonzero(npAND(g2part,behavData.trialsEachBlock[:,3]))
                g2 = np.flatnonzero(npAND(g3part,behavData.trialsEachBlock[:,2]))
                g4 = np.flatnonzero(npAND(g3part,behavData.trialsEachBlock[:,4]))
            else:
                trialsEachCondLabels = ['14kHz (L)','14kHz (R)','14kHz (L)','14kHz (R)']
                colorEachCond = [cp.TangoPalette['ScarletRed1'], cp.TangoPalette['Orange2'],
                                 cp.TangoPalette['ScarletRed1'], cp.TangoPalette['Orange2']]
                #g1part = npAND(npAND(behavData.lowFreqs,behavData.leftReward),npNOT(behavData.early))
                #g2part = npAND(npAND(behavData.lowFreqs,behavData.rightReward),npNOT(behavData.early))
                #g3part = npAND(npAND(behavData.highFreqs,behavData.leftReward),npNOT(behavData.early))
                #g4part = npAND(npAND(behavData.highFreqs,behavData.rightReward),npNOT(behavData.early))
                g1 = np.flatnonzero(npAND(g3part,behavData.trialsEachBlock[:,1]))
                g3 = np.flatnonzero(npAND(g3part,behavData.trialsEachBlock[:,3]))
                g2 = np.flatnonzero(npAND(g2part,behavData.trialsEachBlock[:,2]))
                g4 = np.flatnonzero(npAND(g2part,behavData.trialsEachBlock[:,4]))
        else:
            g1 = np.flatnonzero(npAND(behavData.lowFreqs,behavData.leftReward))
            g2 = np.flatnonzero(npAND(behavData.lowFreqs,behavData.rightReward))
            g3 = np.flatnonzero(npAND(behavData.highFreqs,behavData.leftReward))
            g4 = np.flatnonzero(npAND(behavData.highFreqs,behavData.rightReward))
        # -- Sort if requested --
        if sortby.size:
            # Test for NaN (since numpy<1.4 cannot sort NaNs)
            if np.any(np.isnan(sortby[np.concatenate([g1,g2,g3,g4])])):
                raise ValueError('There are NaNs in sorting vector.')
            g1 = g1[np.argsort(sortby[g1])]
            g2 = g2[np.argsort(sortby[g2])]
            g3 = g3[np.argsort(sortby[g3])]
            g4 = g4[np.argsort(sortby[g4])]

        trialsEachCond = [g1,g2,g3,g4]
        nTrialsEachCond = [len(x) for x in trialsEachCond]
        nCond = len(nTrialsEachCond)
        indsEachCond = [np.arange(nTrialsEachCond[0])]
        lastTrialEachCond = np.cumsum(nTrialsEachCond)
        firstTrialEachCond = np.r_[0,lastTrialEachCond[0:-1]]
        for indc in range(1,nCond):
            theseInds = np.arange(nTrialsEachCond[indc])+lastTrialEachCond[indc-1]
            indsEachCond.append(theseInds)
        #yLabelStr = 'Correct trials (sorted by RT)'+'\n'+'LOWleft | LOWright | HIGHleft | HIGHright'
        yLabelStr = 'Trial'
    elif CONDCASE==2:
        pass
    condInfo = {'nTrialsEachCond':nTrialsEachCond,
                'indsEachCond':indsEachCond,
                'colorEachCond':colorEachCond,
                'trialsEachCondLabels':trialsEachCondLabels,
                'firstTrialEachCond':firstTrialEachCond,
                'lastTrialEachCond':lastTrialEachCond,
                'yLabelStr':yLabelStr}
    return (trialsEachCond,condInfo)


def align_to_event(behavData,CASE):
    if CASE==1:
        xLabelStr = 'Time from sound onset'
        eventOfInterest = behavData.targetOnsetTime - \
                          behavData.trialStartTime + \
                          behavData.trialStartTimeEphys
        '''
        eventOfInterest = behavData.targetOnsetTime[trialsOfInterest] - \
                          behavData.trialStartTime[trialsOfInterest] + \
                          behavData.trialStartTimeEphys[trialsOfInterest]
        '''
    elif CASE==2:
        xLabelStr = 'Time from center poke out'
        eventOfInterest = behavData.centerOutTime - \
                          behavData.trialStartTime + \
                          behavData.trialStartTimeEphys
    elif CASE==3:
        xLabelStr = 'Time from reward port in'
        eventOfInterest = behavData.sideInTime - \
                          behavData.trialStartTime + \
                          behavData.trialStartTimeEphys
        '''
    elif CASE==4:
        xLabelStr = 'UNSORTED Time from side port in'
        eventOfInterest = behavData.sideInTime - \
                          behavData.trialStartTime + \
                          behavData.trialStartTimeEphys
    elif CASE==5:
        xLabelStr = 'UNSORTED Time from center poke out'
        eventOfInterest = behavData.centerOutTime - \
                          behavData.trialStartTime + \
                          behavData.trialStartTimeEphys
    elif CASE==6:
        xLabelStr = 'UNSORTED Time from sound onset (ms)'
        '''
    return (eventOfInterest,xLabelStr)
    
def save_data_each_cell(cellDB,outputDir,timeRange=np.array([-0.3,0.9]),lockTo=1):
    allPostfix = {1:'SoundOn',2:'Cout',3:'SideIn'} ### FIXME: HARDCODED !!!
    for indcell,onecell in enumerate(cellDB):
        #cellStr = '%s_%s_T%dc%d'%(onecell.animalName, onecell.ephysSession,
        #                          onecell.tetrode, onecell.cluster)
        cellStr = str(onecell).replace(' ','_')
        try:
            (behavData,trialEvents,dataTT,spikeInds) = load_cell_reversal(onecell)
        except IOError:
            print 'WARNING: File not found for cell %s'%cellStr
            continue
        (eventOfInterest,xLabelStr) = align_to_event(behavData,lockTo)
        # -- Ignore trialsToExclude --
        eventOfInterest[onecell.trialsToExclude] = np.nan
        ##if len(onecell.trialsToExclude)>0:
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(dataTT.timestamps[spikeInds],
                                                  eventOfInterest,timeRange)
        filePostfix = allPostfix[lockTo]
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        fileName = os.path.join(outputDir,cellStr+'_'+filePostfix+'.npz')
        print fileName
        np.savez(fileName, spikeTimesFromEventOnset=spikeTimesFromEventOnset,
                 trialIndexForEachSpike=trialIndexForEachSpike,
                 indexLimitsEachTrial=indexLimitsEachTrial, xLabelStr=xLabelStr,
                 timeRange=timeRange,animalName=onecell.animalName,
                 behavSession=onecell.behavSession)

def save_data_each_mu(muDB,outputDir,timeRange=np.array([-0.3,0.9]),lockTo=1):
    for indmu,onemu in enumerate(muDB):
        muStr = '%s_%s_T%dmu'%(onemu.animalName, onemu.ephysSession,
                               onemu.tetrode)
        (behavData,trialEvents,dataTT,spikeInds) = load_mu_reversal(onemu)
        (eventOfInterest,xLabelStr) = align_to_event(behavData,lockTo)
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(dataTT.timestamps[spikeInds],
                                                  eventOfInterest,timeRange)
        fileName = os.path.join(outputDir,muStr+'_'+filePostfix+'.npz')
        print fileName
        np.savez(fileName, spikeTimesFromEventOnset=spikeTimesFromEventOnset,
                 trialIndexForEachSpike=trialIndexForEachSpike,
                 indexLimitsEachTrial=indexLimitsEachTrial,
                 timeRange=timeRange,animalName=onemu.animalName,
                 behavSession=onemu.behavSession)

'''
run -e allcells_saja099.py
outputDir='/var/data/neuralynx/saja099_processed'
from extracellpy import sessionanalysis
time sessionanalysis.save_data_each_cell(cellDB,outputDir)

CPU times: user 134.40 s, sys: 36.70 s, total: 171.09 s

run -e allcells_saja099.py
from extracellpy import sessionanalysis
(behavData,trialEvents,dataTT,spikeInds) = sessionanalysis.load_cell_reversal(cellDB[0])

'''


def evaluate_response_each_cell(cellDB,baseRange,rangeStart):
    timeRange = np.array([baseRange[0], baseRange[-1]+np.diff(baseRange)])
    for indcell,onecell in enumerate(cellDB):
        print '[%d/%d] %s'%(indcell+1,nCells,onecell)
        (behavData,trialEvents,dataTT,spikeInds) = sessionanalysis.load_cell_reversal(onecell)
        (eventOfInterest,xLabelStr) = sessionanalysis.align_to_event(behavData,1)
        (trialsEachCond,condInfo) = sessionanalysis.trials_by_condition(behavData,1,outcome='correct')
        trialsOfInterest = np.hstack(trialsEachCond)
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(dataTT.timestamps[spikeInds],
                                                  eventOfInterest[trialsOfInterest],
                                                  timeRange)
        nCond = len(trialsEachCond)
        zStats = np.empty((len(rangeStart),nCond))
        for indCond in range(nCond):
            toUse = np.arange(condInfo['firstTrialEachCond'][indCond],
                              condInfo['lastTrialEachCond'][indCond])
            (zStats[:,indCond],pValues) = spikesanalysis.evaluate_responsiveness(spikeTimesFromEventOnset,
                                                   indexLimitsEachTrial[:,toUse],baseRange,
                                                   rangeStart)
        zStatsEachCell[:,:,indcell] = zStats
        return zStatsEachCell

def evaluate_modulation_each_cell(cellDB,responseRange):
    '''
    '''

    ############# FINISH THIS (see spikesanalysis.evaluate_modulation()  ##########

    timeRange = np.array([baseRange[0], baseRange[-1]+np.diff(baseRange)])
    for indcell,onecell in enumerate(cellDB):
        print '[%d/%d] %s'%(indcell+1,nCells,onecell)
        (behavData,trialEvents,dataTT,spikeInds) = sessionanalysis.load_cell_reversal(onecell)
        (eventOfInterest,xLabelStr) = sessionanalysis.align_to_event(behavData,1)
        (trialsEachCond,condInfo) = sessionanalysis.trials_by_condition(behavData,1,outcome='correct')
        trialsOfInterest = np.hstack(trialsEachCond)
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(dataTT.timestamps[spikeInds],
                                                  eventOfInterest[trialsOfInterest],
                                                  timeRange)
        nCond = len(trialsEachCond)
        zStats = np.empty((len(rangeStart),nCond))
        for indCond in range(nCond):
            toUse = np.arange(condInfo['firstTrialEachCond'][indCond],
                              condInfo['lastTrialEachCond'][indCond])
            (zStats[:,indCond],pValues) = spikesanalysis.evaluate_responsiveness(spikeTimesFromEventOnset,
                                                   indexLimitsEachTrial[:,toUse],baseRange,
                                                   rangeStart)
        zStatsEachCell[:,:,indcell] = zStats
        return zStatsEachCell
'''
    np.savez('/var/tmp/saja099_responsiveness.npz',zStatsEachCell=zStatsEachCell,
         rangeStart=rangeStart,baseRange=baseRange,cellDB=allcells.cellDB)
'''

#def findOnset(boolArray):
#    '''Find indexes where value changes from False to True.'''
############ SEE INSTEAD METHOD IN LoadNeuralynx.EVENTS #######

def find_photostim_times(oneCell,bitTRIALIND=bitTRIALIND_DEFAULT,npulses=1):
    '''
    Find times between trialStart and photoStim (from ephys data)
    '''
    animalName = oneCell.animalName
    ephysSession = oneCell.ephysSession
    dataDir = os.path.join(EPHYSPATH,'%s/%s/'%(animalName,ephysSession))
    eventsFile = os.path.join(dataDir,'Events.nev')
    events = loadneuralynx.DataEvents(eventsFile)
    trialEvents = (events.valueTTL & (1<<bitTRIALIND)) != 0
    trialStartTimeNL = 1e-6*events.timestamps[trialEvents]

    (trialOnset,trialOffset) = events.find_bit_changes(bitTRIALIND)
    (targetOnset,targetOffset) = events.find_bit_changes(bitTARGETIND)
    (photoOnset,photoOffset) = events.find_bit_changes(bitPHOTOSTIMIND)

    eventsTime={}
    eventsTime['trial'] = 1e-6*events.timestamps[trialOnset]
    eventsTime['target'] = 1e-6*events.timestamps[targetOnset]
    eventsTime['photostim'] = 1e-6*events.timestamps[photoOnset]


    ######### VERIFY THIS ############
    if npulses>1:
        eventsTime['photostim'] = eventsTime['photostim'][0::npulses]
    # -- Exclude second photostim --
    #print('WARNING: using all photostim events (not correct if trains were presented)')
    #eventsTime['photostim'] = eventsTime['photostim'][0::2]

    ########## DEBUG ###########
    #eventsTime['photostim'] = eventsTime['target']
    
    # --- Find time of photostim for each trial ---
    nTrialsEphys = len(eventsTime['trial'])
    nTarget = len(eventsTime['target'])
    nPhotoStim = len(eventsTime['photostim'])
    diffTargetPhoto = np.tile(eventsTime['photostim'],(nTrialsEphys,1)).T - \
                      np.tile(eventsTime['trial'],(nPhotoStim,1))
    diffTargetPhoto[diffTargetPhoto<0]=1e6 # Arbitrary large number
    trialForEachPhoto = np.argmin(diffTargetPhoto,axis=1)

    photoStimOnset = np.empty(nTrialsEphys)
    photoStimOnset.fill(np.nan)
    photoStimOnset[trialForEachPhoto] = eventsTime['photostim']-eventsTime['trial'][trialForEachPhoto]
    # Add empty trial at the beginning (to account for BControl empty trial)
    photoStimOnset = np.hstack((np.nan,photoStimOnset[:-1]))
    return photoStimOnset


if __name__ == "__main__":
    import celldatabase
    reload(celldatabase)

    CASE = 3

    if CASE == 1:
        import pylab as plt
        animalName = 'saja099'
        ephysSession = '2011-04-18_20-44-49'
        behavSession = '20110418a'
        tetrode = 1
        cluster = 3
        onecell = celldatabase.CellInfo(animalName,ephysSession,behavSession,tetrode,cluster)
        print 'Loading %s ...'%(onecell)
        (behavData,trialEvents,dataTT,spikeInds) = load_cell_reversal(onecell)
        timeRange = np.array([-0.3,0.9])
        trialsOfInterest = range(300)
        (eventOfInterest,xLabelStr) = align_to_event(behavData,1)
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(dataTT.timestamps[spikeInds],
                                                  eventOfInterest,timeRange)
        pR, = plt.plot(1e3*spikeTimesFromEventOnset,trialIndexForEachSpike,'.k')
        pR.set_markersize(2)
        plt.draw()
        plt.show()
    elif CASE==2:
        animalName   = 'saja100'
        ephysSession = '2011-10-26_11-57-31'
        behavSession = '20111026a'
        tetrode = 1
        clusters = [2,3,4,5,6,7,9,10,11]
        oneSite = celldatabase.MultiUnitInfo(animalName = animalName,
                                             ephysSession = ephysSession,
                                             behavSession = behavSession,
                                             tetrode = tetrode,
                                             clusters = clusters)
        (behavData,trialEvents,dataTT,spikeInds) = load_mu_reversal(oneSite)
    elif CASE==3:
        '''Test loading LFP with behavior '''
        oneLFP = celldatabase.LFPInfo('saja125','2012-02-02_16-33-29','20120202a',17)
        (behavData,trialEvents,dataLFP) = load_lfp_reversal(oneLFP)
