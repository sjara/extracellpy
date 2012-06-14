#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Functions for analysis of spikes
'''

#from pylab import *
import numpy as np
import os

__author__ = 'Santiago Jaramillo'
__version__ = '0.2'

npAND = np.logical_and
npOR = np.logical_or

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


def eventlocked_spiketimes(timeStamps,eventOnsetTimes,timeRange,spikeindex=False):
    '''Create a vector with the spike timestamps w.r.t. event onset.

    (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = 
        eventlocked_spiketimes(timeStamps,eventOnsetTimes,timeRange)

    timeStamps: the time of each spike.
    eventOnsetTimes: the time of each instance of the event to lock to.
    timeRange: two-element array specifying time-range to extract around event.

    SpikeTimesFromEventOnset: 1D array with time of spikes locked to event.
    TrialIndexForEachSpike: 1D array with the trial corresponding to each spike.
    IndexLimitsEachTrial: [2,nTrials] range of spikes for each trial. Note that
       the range is from firstSpike to lastSpike+1 (like in python slices)

    TrialIndexForEachSpike starts at 0 (not 1 like matlab)
    '''
    nTrials = len(eventOnsetTimes)
    # FIXME: check if timestamps are sorted
    
    SpikeTimesFromEventOnset = np.empty(0,dtype='float64')
    spikeIndices = np.empty(0,dtype='int')
    TrialIndexForEachSpike = np.empty(0,dtype='int')
    IndexLimitsEachTrial = np.empty((2,nTrials),dtype='int')
    accumIndexFirstSpike = 0
    
    for indtrial in np.arange(nTrials):
        ThisTrialRange = eventOnsetTimes[indtrial] + timeRange
        FirstSpikeInTrial = np.searchsorted(timeStamps,ThisTrialRange[0])
        #LastSpikeInTrial = np.searchsorted(timeStamps,ThisTrialRange[-1])-1
        ####### FIX ME: LastSpikeInTrial should not be negative, because
        ### slice(0,-1) means from 0 to lastIndex
        #SpikesThisTrial = np.arange(FirstSpikeInTrial,LastSpikeInTrial+1)

        LastSpikeInTrialPlusOne = np.searchsorted(timeStamps,ThisTrialRange[-1])
        #SpikesThisTrial = slice(FirstSpikeInTrial,LastSpikeInTrialPlusOne) # FIXME: Faster?
        SpikesThisTrial = np.arange(FirstSpikeInTrial,LastSpikeInTrialPlusOne)
        NSpikesThisTrial = LastSpikeInTrialPlusOne - FirstSpikeInTrial

        spikeIndices = np.concatenate((spikeIndices,SpikesThisTrial))
        SpikeTimesFromEventOnset = np.concatenate((SpikeTimesFromEventOnset,
                                        timeStamps[SpikesThisTrial]-eventOnsetTimes[indtrial]))
        TrialIndexForEachSpike = np.concatenate((TrialIndexForEachSpike,
                                            np.repeat(indtrial,NSpikesThisTrial)))
        #IndexLimitsEachTrial[:,indtrial] = [accumIndexFirstSpike,accumIndexFirstSpike+NSpikesThisTrial-1] # OLD:non python-slice
        IndexLimitsEachTrial[:,indtrial] = [accumIndexFirstSpike,accumIndexFirstSpike+NSpikesThisTrial]
        accumIndexFirstSpike += NSpikesThisTrial
        #1/0 ### DEBUG
    if spikeindex:
        return (SpikeTimesFromEventOnset,TrialIndexForEachSpike,IndexLimitsEachTrial,spikeIndices)
    else:
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

#@profile
def calculate_psth(spikeRasterMat,timeVec,windowSize):
    '''Calculate Peri-Stimulus Time Histogram.
    It uses a square window and returns the average spikes per second.
    '''
    nTrials = spikeRasterMat.shape[0]
    deltaTime = timeVec[1]-timeVec[0]
    windowSizeInSamples = int(round(windowSize/deltaTime))

    windowShape = np.ones(windowSizeInSamples,dtype=np.float64)
    windowShape = windowShape/(windowSizeInSamples*deltaTime)
    '''
    PSTHeach = np.empty(spikeRasterMat.shape,dtype=np.float64)
    for indt,trial in enumerate(spikeRasterMat):
        PSTHeach[indt,:] = np.convolve(trial,windowShape,'same')
    PSTH = np.mean(PSTHeach,axis=0)
    '''
    spikeMatMean = np.mean(spikeRasterMat,axis=0)
    PSTH = np.convolve(spikeMatMean,windowShape,'same')
    return PSTH

def calculate_psth_per_condition(spikeTimesFromEventOnset,indexLimitsEachTrial,
                                 timeVec,windowSizePSTH, trialsEachCond):
    '''Calculate PSTH for each condition defined by trialsEachCond.
    See: 'spiketimes_to_sparsemat' and 'calculate_psth' for details
         on how PSTH is calculated.
    '''
    spikeRasterMat = spiketimes_to_sparsemat(spikeTimesFromEventOnset,
                                             indexLimitsEachTrial,
                                             timeVec)
    PSTH = np.empty((len(trialsEachCond),spikeRasterMat.shape[1]))
    for (indc,indsThisCond) in enumerate(trialsEachCond):
        thisSpikeMat = spikeRasterMat[indsThisCond,:]
        PSTH[indc,:] = calculate_psth(thisSpikeMat, timeVec, windowSizePSTH)
    binsStartTime = timeVec
    return(PSTH,binsStartTime,spikeRasterMat)


def spiketimes_to_sparsemat(spikeTimesFromEventOnset,indexLimitsEachTrial,timeVec):
    '''Create a (sparse) matrix with spikes given the times of the spikes.

       spikeTimesFromEventOnset: vector of spikes timestamps with respect
         to the onset of the event.
       indexLimitsEachTrial: each column contains [firstInd,lastInd+1] of the spikes on a trial
       timeVec: vector that defines the horizontal axis of the matrix.

       Returns spikeRasterMat
       NOTES: 
       - If two spikes are within the same bin, they are counted as only one.
    '''
    nTrials = indexLimitsEachTrial.shape[1]
    deltaTime = timeVec[1]-timeVec[0]
    # NOTE: Using full matrix instead of sparse.
    spikeRasterMat = np.empty((nTrials,len(timeVec)),dtype=bool)
    for indtrial in range(nTrials):
        indsThisTrial = slice(indexLimitsEachTrial[0,indtrial],indexLimitsEachTrial[1,indtrial])
        sampleIndexOfSpike = np.around( (spikeTimesFromEventOnset[indsThisTrial]-timeVec[0]) / deltaTime ).astype(int)
        #if indtrial==1: 1/0 ### DEBUG
        spikeRasterMat[indtrial,sampleIndexOfSpike]=True
    return spikeRasterMat


def count_spikes_in_range(spikeTimesFromEventOnset,indexLimitsEachTrial,timeRange):
    '''Count number of spikes on each trial in a given time range.

       spikeTimesFromEventOnset: vector of spikes timestamps with respect
         to the onset of the event.
       indexLimitsEachTrial: each column contains [firstInd,lastInd+1] of the spikes on a trial
       timeRange: time range to evaluate

       returns nSpikes
    '''
    nTrials = indexLimitsEachTrial.shape[1]
    nSpikes = np.empty(nTrials)
    for indtrial in range(nTrials):
        indsThisTrial = slice(indexLimitsEachTrial[0,indtrial],indexLimitsEachTrial[1,indtrial])
        #sampleIndexOfSpike = np.around( (spikeTimesFromEventOnset[indsThisTrial]-timeVec[0]) / deltaTime ).astype(int)
        spikeTimesThisTrial = spikeTimesFromEventOnset[indsThisTrial]
        nSpikes[indtrial] = sum(npAND(spikeTimesThisTrial>timeRange[0],spikeTimesThisTrial<timeRange[-1]))
    return nSpikes
    

def OLD_spiketimes_to_sparsemat(spikeTimesFromEventOnset,trialIndexForEachSpike,
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

def evaluate_responsiveness(spikeTimesFromEventOnset,indexLimitsEachTrial,baseRange,
                            rangeStart):
    '''Evaluate the probability of observing firing for given periods given the
       firing of a base range.
    '''
    from scipy import stats
    rangeLength = np.diff(baseRange)
    pValueEachRange = np.empty(len(rangeStart))
    zStatsEachRange = np.empty(len(rangeStart))
    nspkBase=count_spikes_in_range(spikeTimesFromEventOnset,
                                   indexLimitsEachTrial,baseRange)
    for indbin,thisStart in enumerate(rangeStart):
        respRange = np.r_[0,rangeLength] + thisStart
        nspkResp=count_spikes_in_range(spikeTimesFromEventOnset,
                                       indexLimitsEachTrial,respRange)
        #[tStat,pValue] = stats.ttest_ind(nspkResp,nspkBase)
        [zStat,pValue] = stats.ranksums(nspkResp,nspkBase)
        pValueEachRange[indbin] = pValue
        zStatsEachRange[indbin] = zStat
    return (zStatsEachRange,pValueEachRange)


def evaluate_modulation(spikeTimesFromEventOnset,indexLimitsEachTrial,responseRange,trialsEachCond):
    '''Evaluate the response for each of two conditions and the probability of observing this
       difference by chance.

       spikeTimesFromEventOnset
       indexLimitsEachTrial
       responseRange
       trialsEachCond: list of two arrays of indexes (either bool or int)

       Returns:
       (meanSpikes,pValue)
    '''
    from scipy import stats
    nspkResp = []
    for indcond,trialsThisCond in enumerate(trialsEachCond):
        indexLimitsSubset = indexLimitsEachTrial[:,trialsThisCond]
        nspkResp.append(count_spikes_in_range(spikeTimesFromEventOnset,
                                              indexLimitsSubset,responseRange))
    meanSpikes = np.array([np.mean(n) for n in nspkResp])
    [zStat,pValue] = stats.ranksums(nspkResp[0],nspkResp[1])
    return (meanSpikes,pValue)


def estimate_spike_shape(onecell):
    '''
    Estimate average waveform of spike and some measurements.
    '''
    import settings
    import loadneuralynx
    
    animalName = onecell.animalName
    ephysSession = onecell.ephysSession
    behavSession = onecell.behavSession
    tetrode = onecell.tetrode
    cluster = onecell.cluster
    cellStr = str(onecell).replace(' ','_')

    dataDir = os.path.join(settings.EPHYS_PATH,'%s/%s/'%(animalName,ephysSession))
    clustersDir = os.path.join(settings.EPHYS_PATH,'%s/%s_kk/'%(animalName,ephysSession))

    tetrodeFile = os.path.join(dataDir,'TT%d.ntt'%tetrode)
    #if not prevCell or onecell.tetrode!=prevCell.tetrode:
    dataTT = loadneuralynx.DataTetrode(tetrodeFile,readWaves=True)
    newShape = (dataTT.params['NumADChannels'],dataTT.params['WaveformLength'],-1)
    dataTT.samples = dataTT.samples.reshape(newShape,order='F')
    # -- Load clusters --
    clusterFileSuffix = '1'
    clustersFile = os.path.join(clustersDir,'TT%d.clu.%s'%(tetrode,clusterFileSuffix))
    dataTT.set_clusters(clustersFile)
    spikeInds = np.flatnonzero(dataTT.clusters==cluster)

    nWavesToAverage = min(len(spikeInds),2000)
    peakSample = np.argmax(dataTT.samples[:,:,:nWavesToAverage],axis=1)
    commonPeakSample = np.bincount(peakSample.flatten()).argmax()
    avWaveforms = np.mean(dataTT.samples[:,:,:nWavesToAverage],axis=2)
    maxChannel = np.argmax(avWaveforms[:,commonPeakSample])
    selectedWaveforms = peakSample[maxChannel,:]==commonPeakSample
    waveform = np.mean(dataTT.samples[maxChannel,:,np.flatnonzero(selectedWaveforms)],axis=0)
    # NOTE: indexing seems to change the order or things, so we need to take the mean across axis=0.
    spikeMeasures = {}
    spikeMeasures['peakSample'] = np.argmax(waveform)
    spikeMeasures['valleySample'] = np.argmin(waveform)
    spikeMeasures['spikeWidth'] = spikeMeasures['valleySample']-spikeMeasures['peakSample']
    spikeMeasures['maxValue'] = np.max(waveform)
    spikeMeasures['minValue'] = np.min(waveform)
    #plot(waveform)
    return (waveform,spikeMeasures)




def load_mclust_t(fileName):
    from struct import unpack
    fid = open(fileName,'rb')
    # Read header
    while(fid.readline()!='%%ENDHEADER\n'): pass
    packedRecords = fid.read()
    fid.close()
    nRec = len(packedRecords)/4 # Because they are uint32
    timeStamps = unpack('>'+'L'*nRec,packedRecords)
    return 1e-4*np.array(timeStamps)


'''
# THE PROBLEM WITH THIS IS THE DEPENDENCY ON load_mu_reversal()
# IT IS NOT GENERAL TO ANY EPHYS DATA
def save_locked_spikes():
    (behavData,trialEvents,dataTT,spikeInds) = load_mu_reversal(onemu)
    (eventOfInterest,xLabelStr) = align_to_event(behavData,1)
    (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(dataTT.timestamps[spikeInds],
                                              eventOfInterest,timeRange)
    fileName = os.path.join(outputDir,muStr+'.npz')
    print fileName
    np.savez(fileName, spikeTimesFromEventOnset=spikeTimesFromEventOnset,
             trialIndexForEachSpike=trialIndexForEachSpike,
             indexLimitsEachTrial=indexLimitsEachTrial,
             timeRange=timeRange,animalName=onemu.animalName,
             behavSession=onemu.behavSession)
'''

#if __name__ == "__main__":

'''
=============== TEST CODE ==============
    
# -- Test eventlocked_spiketimes --
from extracellpy import spikesanalysis
timeStamps = np.array([0.1,0.2,0.3, 1.1,1.3, 2.2, 3.5, 5.1])
eventOnsetTimes = np.array([-2,-1,0,1,2,3,4,5])
timeRange = [-0.4,0.4]
(SpikeTimesFromEventOnset,TrialIndexForEachSpike,IndexLimitsEachTrial) = spikesanalysis.eventlocked_spiketimes(timeStamps,eventOnsetTimes,timeRange)
print np.c_[SpikeTimesFromEventOnset,TrialIndexForEachSpike]
print IndexLimitsEachTrial


'''
