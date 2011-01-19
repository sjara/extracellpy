#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Functions for analysis of spikes
'''

from pylab import *

__author__ = 'Santiago Jaramillo'
__version__ = "$Revision: 2009-03-06$"


def findfirstsmallerthan(xvec,thevalue):
    '''Find index for first value in xvec smaller than thevalue
    '''
    ind=None
    fromvalue=0
    tovalue=len(xvec)
    while (tovalue-fromvalue)>0:
        mid = ceil((fromvalue + tovalue)/2.0)
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


def eventlocked_spiketimes(TimeStamps,EventOnsetTimes,TimeRange):
    '''Create a vector with the spike timestamps w.r.t. event onset.

       TrialIndexForEachSpike starts at 0 (not 1 like matlab)
    '''
    #[SpikeTimesFromEventOnset,TrialIndexForEachSpike] = ...

    nTrials = len(EventOnsetTimes)

    ### CHECK IF TIMESTAMPS ARE SORTED ###
    
    SpikeTimesFromEventOnset = []
    TrialIndexForEachSpike = []

    for indtrial in arange(nTrials):
        ThisTrialRange = EventOnsetTimes[indtrial] + TimeRange
        FirstSpikeInTrial = findfirstsmallerthan(TimeStamps,ThisTrialRange[0])+1
        LastSpikeInTrial = findfirstsmallerthan(TimeStamps,ThisTrialRange[-1])
        SpikesThisTrial = arange(FirstSpikeInTrial,LastSpikeInTrial+1)
        SpikeTimesFromEventOnset = r_[SpikeTimesFromEventOnset,\
                         TimeStamps[SpikesThisTrial]-EventOnsetTimes[indtrial]]
        NSpikesThisTrial = len(SpikesThisTrial)
        TrialIndexForEachSpike = r_[TrialIndexForEachSpike,\
                                    repeat(indtrial,NSpikesThisTrial)];
    
    return (SpikeTimesFromEventOnset,TrialIndexForEachSpike)


def calculate_psth_quant(SpikeRasterMat,TimeVec,WindowSize,TimeOffset=0):
    '''Calculate Peri-Stimulus Time Histogram without overlapping windows.
    '''
    Ntrials = SpikeRasterMat.shape[0]
    deltaTime = TimeVec[1]-TimeVec[0]
    WindowSizeInSamples = int(round(WindowSize/deltaTime))
    OffsetInSamples = int(round(TimeOffset/deltaTime))

    BinsStartSample = arange(OffsetInSamples,len(TimeVec)-OffsetInSamples,WindowSizeInSamples)
    BinsStartTime = BinsStartSample*deltaTime + TimeVec[0]

    WindowCount = zeros((Ntrials,len(BinsStartSample)))
    for indb in range(len(BinsStartSample)-1):
        WindowCount[:,indb] = sum(SpikeRasterMat[:,\
                                  arange(BinsStartSample[indb],BinsStartSample[indb+1])] ,axis=1)
    PSTH = mean(WindowCount,axis=0)
    deltaBinTime = BinsStartTime[1]-BinsStartTime[0]
    PSTH = PSTH/deltaBinTime
    #return (PSTH,BinsStartTime,WindowCount)
    return (PSTH,BinsStartTime)


def spiketimes_to_sparsemat(SpikeTimesFromEventOnset,TrialIndexForEachSpike,\
                                         Ntrials, TimeVec):
    '''Create a (sparse) matrix with spikes given the times of the spikes.

       SpikeTimesFromEventOnset: vector of spikes timestamps with respect
       to the onset of the event.
       TrialIndexForEachSpike: trial index for each spike (starting at 0, not 1 like matlab)
       NTrials: total number of trials (necessary since some trials may have no spikes)
       TimeVec: vector that defines the horizontal axis of the matrix.

       Returns SpikeRasterMat
       NOTES: 
       - If two spikes are within the same bin, they are counted as only one.

    Based on spiketimes_to_sparsemat.m (by Santiago Jaramillo - 2007.05.09)
    '''
    deltaTime = TimeVec[1]-TimeVec[0]
    SpikeRasterMat = zeros((Ntrials,len(TimeVec)),dtype=bool)
    #print('Using full matrix instead of sparse.')
    for indtrial in range(Ntrials):
        SpikeTimesThisTrial = SpikeTimesFromEventOnset[TrialIndexForEachSpike==indtrial]
        SampleIndexOfSpike = around( (SpikeTimesThisTrial-TimeVec[0]) / deltaTime ).astype(int)
        SpikeRasterMat[indtrial,SampleIndexOfSpike]=True

    return SpikeRasterMat

def sparsemat_to_spiketimes(SpikeRasterMat, TimeVec):
    '''Calculate the spike times given a (sparse) matrix and time vector.
       It should be used only for displaying data since it returns
       quantized times.

       SpikeTimesFromEventOnset: vector of spikes timestamps with respect
         to the onset of the event.
       TrialIndexForEachSpike: trial index for each spike.
       TimeVec: vector that defines the horizontal axis of the matrix.

       Returns (SpikeTimesFromEventOnset,TrialIndexForEachSpike)
       
       Based on sparsemat_to_spiketimes.m (by Santiago Jaramillo - 2007.08.09)
    '''
    (TrialIndexForEachSpike,Cols) = SpikeRasterMat.nonzero()
    SpikeTimesFromEventOnset = TimeVec[Cols]
    return (SpikeTimesFromEventOnset,TrialIndexForEachSpike)
