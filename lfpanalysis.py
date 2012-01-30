#!/usr/bin/env python
# -*- coding: utf-8 -*-

############# OBSOLETE #############
### see loadneuralynx.DataCont instead ###

'''Functions for analysis of local field potentials (LFP)
'''

#from pylab import *
import numpy as np

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

npAND = np.logical_and
npOR = np.logical_or

def eventlocked_lfp(dataLFP,eventOnsetTimes,timeRange):
    '''Make matrix of LFP traces locked to stimulus
    dataLFP should be created by loadneuralynx and have
    dataLFP.samplingRate, dataLFP.samples, dataLFP.time
    '''
    wintoplot=np.arange(-20000,40000)
    timeVec = np.arange(timeRange[0],1/dataLFP.samplingRate,timeRange[-1])
    firstSample = round(timeRange[0]*dataLFP.samplingRate)
    nSamples = len(timeVec)
    nTrials = len(eventOnsetTimes)
    samplesRange = np.arange(nSamples)+firstSample
    lockedLFP = np.empty((nTrials,nSamples))
    for inde,eventTime in enumerate(eventOnsetTimes):
        eventSample = np.searchsorted(dataLFP.time,eventTime)
        thiswin = samplesRange+eventSample
        lockedLFP[inde,:] = dataLFP.samples[thiswin]
    return (lockedLFP,timeVec)
