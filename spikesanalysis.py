#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Functions for analysis of spikes
'''

#from pylab import *
import numpy as np

__author__ = 'Santiago Jaramillo'
__version__ = '0.2'


def findfirstsmallerthan(xvec,thevalue):
    '''Find index for first value in xvec smaller than thevalue
       OBSOLETE: use np.searchsorted instead
    '''
    ind=None
    fromvalue=0
    tovalue=len(xvec)
    while (tovalue-fromvalue)>0:
        mid = np.ceil((fromvalue + tovalue)/2.0)
        #mid = (fromvalue + tovalue)/2
        if(xvec[mid]>thevalue):
            tovalue = mid-1
        else:
            fromvalue = mid
    # -- Check if first value in vector is smaller than thevalue --
    if xvec[tovalue]<thevalue:
        ind=tovalue.astype(np.int)
    else:
        ind=None
        #raise ValueError('SOME ERROR HERE')
    return ind


def OLD_eventlocked_spiketimes(TimeStamps,EventOnsetTimes,TimeRange):
    '''Create a vector with the spike timestamps w.r.t. event onset.

    [SpikeTimesFromEventOnset,TrialIndexForEachSpike] = 
        eventlocked_spiketimes(TimeStamps,EventOnsetTimes,TimeRange)

    TrialIndexForEachSpike starts at 0 (not 1 like matlab)
    '''
    nTrials = len(EventOnsetTimes)
    # FIXME: check if timestamps are sorted
    
    SpikeTimesFromEventOnset = np.empty(0,dtype='float64')
    TrialIndexForEachSpike = np.empty(0,dtype='int')

    for indtrial in np.arange(nTrials):
        ThisTrialRange = EventOnsetTimes[indtrial] + TimeRange
        FirstSpikeInTrial = np.searchsorted(TimeStamps,ThisTrialRange[0])
        LastSpikeInTrial = np.searchsorted(TimeStamps,ThisTrialRange[-1])-1
        SpikesThisTrial = np.arange(FirstSpikeInTrial,LastSpikeInTrial+1)
        SpikeTimesFromEventOnset = np.concatenate((SpikeTimesFromEventOnset,
                                        TimeStamps[SpikesThisTrial]-EventOnsetTimes[indtrial]))
        NSpikesThisTrial = len(SpikesThisTrial)
        TrialIndexForEachSpike = np.concatenate((TrialIndexForEachSpike,
                                            np.repeat(indtrial,NSpikesThisTrial)))
        #SpikeTimesFromEventOnset = r_[SpikeTimesFromEventOnset,\
        #                 TimeStamps[SpikesThisTrial]-EventOnsetTimes[indtrial]]
        #TrialIndexForEachSpike = r_[TrialIndexForEachSpike,\
        #                            repeat(indtrial,NSpikesThisTrial)];
    
    return (SpikeTimesFromEventOnset,TrialIndexForEachSpike)


def eventlocked_spiketimes(TimeStamps,EventOnsetTimes,TimeRange):
    '''Create a vector with the spike timestamps w.r.t. event onset.

    [SpikeTimesFromEventOnset,IndexLimitsEachTrial] = 
        eventlocked_spiketimes(TimeStamps,EventOnsetTimes,TimeRange)

    TrialIndexForEachSpike starts at 0 (not 1 like matlab)
    '''
    nTrials = len(EventOnsetTimes)
    # FIXME: check if timestamps are sorted
    
    SpikeTimesFromEventOnset = np.empty(0,dtype='float64')
    TrialIndexForEachSpike = np.empty(0,dtype='int')
    IndexLimitsEachTrial = np.empty((2,nTrials),dtype='int')
    accumIndexFirstSpike = 0

    for indtrial in np.arange(nTrials):
        ThisTrialRange = EventOnsetTimes[indtrial] + TimeRange
        FirstSpikeInTrial = np.searchsorted(TimeStamps,ThisTrialRange[0])
        LastSpikeInTrial = np.searchsorted(TimeStamps,ThisTrialRange[-1])-1
        SpikesThisTrial = np.arange(FirstSpikeInTrial,LastSpikeInTrial+1)
        SpikeTimesFromEventOnset = np.concatenate((SpikeTimesFromEventOnset,
                                        TimeStamps[SpikesThisTrial]-EventOnsetTimes[indtrial]))
        NSpikesThisTrial = len(SpikesThisTrial)
        TrialIndexForEachSpike = np.concatenate((TrialIndexForEachSpike,
                                            np.repeat(indtrial,NSpikesThisTrial)))
        IndexLimitsEachTrial[:,indtrial] = [accumIndexFirstSpike,accumIndexFirstSpike+NSpikesThisTrial-1]
        accumIndexFirstSpike += NSpikesThisTrial
        #1/0 ### DEBUG
    return (SpikeTimesFromEventOnset,TrialIndexForEachSpike,IndexLimitsEachTrial)


def calculate_psth_quant(spikeRasterMat,timeVec,windowSize,TimeOffset=0):
    '''Calculate Peri-Stimulus Time Histogram without overlapping windows.'''
    nTrials = spikeRasterMat.shape[0]
    deltaTime = timeVec[1]-timeVec[0]
    windowSizeInSamples = int(round(windowSize/deltaTime))
    OffsetInSamples = int(round(TimeOffset/deltaTime))

    BinsStartSample = np.arange(OffsetInSamples,len(timeVec)-OffsetInSamples,windowSizeInSamples)
    BinsStartTime = BinsStartSample*deltaTime + timeVec[0]

    WindowCount = np.zeros((nTrials,len(BinsStartSample)))
    for indb in range(len(BinsStartSample)-1):
        spikeMatThisBin = spikeRasterMat[:,np.arange(BinsStartSample[indb],BinsStartSample[indb+1])]
        WindowCount[:,indb] = np.sum(spikeMatThisBin,axis=1)
    PSTH = np.mean(WindowCount,axis=0)
    deltaBinTime = BinsStartTime[1]-BinsStartTime[0]
    PSTH = PSTH/deltaBinTime
    #return (PSTH,BinsStartTime,WindowCount)
    return (PSTH,BinsStartTime)

def calculate_psth(spikeRasterMat,timeVec,windowSize):
    ''' '''
    nTrials = spikeRasterMat.shape[0]
    deltaTime = timeVec[1]-timeVec[0]
    windowSizeInSamples = int(round(windowSize/deltaTime))

    windowShape = np.ones(windowSizeInSamples)
    PSTHeach = np.empty(spikeRasterMat.shape,dtype=np.float64)
    for indt,trial in enumerate(spikeRasterMat):
        PSTHeach[indt,:] = np.convolve(trial,windowShape,'same')
    PSTH = np.mean(PSTHeach,axis=0)
    return PSTH


def spiketimes_to_sparsemat(spikeTimesFromEventOnset,indexLimitsEachTrial,\
                                         nTrials, timeVec):
    '''Create a (sparse) matrix with spikes given the times of the spikes.

       spikeTimesFromEventOnset: vector of spikes timestamps with respect
         to the onset of the event.
       indexLimitsEachTrial: each column contains [firstInd,lastInd] of the spikes on a trial
       nTrials: total number of trials (necessary since some trials may have no spikes)
       timeVec: vector that defines the horizontal axis of the matrix.

       Returns spikeRasterMat
       NOTES: 
       - If two spikes are within the same bin, they are counted as only one.
    '''
    deltaTime = timeVec[1]-timeVec[0]
    # NOTE: Using full matrix instead of sparse.
    spikeRasterMat = np.empty((nTrials,len(timeVec)),dtype=bool)
    for indtrial in range(nTrials):
        indsThisTrial = slice(indexLimitsEachTrial[0,indtrial],indexLimitsEachTrial[1,indtrial]+1)
        sampleIndexOfSpike = np.around( (spikeTimesFromEventOnset[indsThisTrial]-timeVec[0]) / deltaTime ).astype(int)
        #if indtrial==1: 1/0 ### DEBUG
        spikeRasterMat[indtrial,sampleIndexOfSpike]=True
    return spikeRasterMat

def OLD_spiketimes_to_sparsemat(spikeTimesFromEventOnset,trialIndexForEachSpike,\
                                         nTrials, timeVec):
    '''Create a (sparse) matrix with spikes given the times of the spikes.

       spikeTimesFromEventOnset: vector of spikes timestamps with respect
       to the onset of the event.
       trialIndexForEachSpike: trial index for each spike (starting at 0, not 1 like matlab)
       NTrials: total number of trials (necessary since some trials may have no spikes)
       timeVec: vector that defines the horizontal axis of the matrix.

       Returns spikeRasterMat
       NOTES: 
       - If two spikes are within the same bin, they are counted as only one.

    Based on spiketimes_to_sparsemat.m (by Santiago Jaramillo - 2007.05.09)
    '''
    deltaTime = timeVec[1]-timeVec[0]
    # NOTE: Using full matrix instead of sparse.
    spikeRasterMat = np.empty((nTrials,len(timeVec)),dtype=bool)
    for indtrial in range(nTrials):
        spikeTimesThisTrial = spikeTimesFromEventOnset[trialIndexForEachSpike==indtrial]
        sampleIndexOfSpike = np.around( (spikeTimesThisTrial-timeVec[0]) / deltaTime ).astype(int)
        #if indtrial==1: 1/0 ### DEBUG
        spikeRasterMat[indtrial,sampleIndexOfSpike]=True
    return spikeRasterMat

def sparsemat_to_spiketimes(spikeRasterMat, timeVec):
    '''Calculate the spike times given a (sparse) matrix and time vector.
       It should be used only for displaying data since it returns
       quantized times.

       spikeTimesFromEventOnset: vector of spikes timestamps with respect
         to the onset of the event.
       trialIndexForEachSpike: trial index for each spike.
       timeVec: vector that defines the horizontal axis of the matrix.

       Returns (spikeTimesFromEventOnset,trialIndexForEachSpike)
       
       Based on sparsemat_to_spiketimes.m (by Santiago Jaramillo - 2007.08.09)
    '''
    (trialIndexForEachSpike,Cols) = spikeRasterMat.nonzero()
    spikeTimesFromEventOnset = timeVec[Cols]
    return (spikeTimesFromEventOnset,trialIndexForEachSpike)
