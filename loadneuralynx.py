#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Objects for reading Neuralynx data.
It depends on Cython implementations from the cLoadNeuralynx module.
'''

import numpy as np
import re
import cLoadNeuralynx

__author__ = 'Santiago Jaramillo'
__version__ = '0.2'

from struct import unpack

class DataCont(cLoadNeuralynx.DataCont):
    '''Access Neuralynx NCS files containing continuous data.
    Reading the data is implemented by the parent class from cLoadNeuralynx.'''
    def __init__(self,fileName):
        cLoadNeuralynx.DataCont.__init__(self,fileName,readTS=True)
        self._test_integrity()
        patt = re.compile(r'-SamplingFrequency (\d+\.\d+)')
        self.samplingRate = float(patt.search(self.header).groups()[0])
        '''
        if len(self.timestamps):
            self.firstTimestamp = self.timestamps[0]
        else:
            self.firstTimestamp = []
        '''
        ###self.time = self._generate_timevec()
    def _test_integrity(self):
        if np.any(np.diff(np.diff(self.timestamps))):
            #raise TypeError('Not all records are contiguous. Packets lost?')
            print('WARNING: Not all records are contiguous. Packets lost?')
    def _generate_timevec(self):
        '''This may be more complicated than necessary: I'm creating an equally
        spaced vector between each record timestamp.
        '''
        dRec = self.timestamps[1]-self.timestamps[0]
        samplesPerRecord = len(self.samples)/self.nRecords
        timeArray = np.tile(self.timestamps.reshape(self.timestamps.shape[0], -1),
                            (1,samplesPerRecord))
        interimTimestamps = np.linspace(0,dRec,samplesPerRecord,endpoint=False).astype(timeArray.dtype)
        timeVector = timeArray+interimTimestamps
        return timeVector.flatten()
    def _parse_header(self):
        pass
    #def __repr__(self):
    #    pass
        '''
        objStrings = []
        for key,value in sorted(vars(self).iteritems()):
            objStrings.append('%s: %s\n'%(key,str(value)))
        return ''.join(objStrings)
        '''
    def time_of_sample(self,sampleInd):
        samplingRateInMicroSec = 1e-6*self.samplingRate
        samplesPerRecord = self.nSamples/self.nRecords
        recordThisSample = sampleInd//samplesPerRecord
        timeThisRecord = self.timestamps[recordThisSample]
        firstSampleThisRecord = samplesPerRecord*recordThisSample
        timeFromRecordStart = (sampleInd-firstSampleThisRecord)/samplingRateInMicroSec
        timeThisSample = timeFromRecordStart + timeThisRecord
        timeThisSample = np.array(timeThisSample,dtype=self.timestamps.dtype)
        return timeThisSample

    def samples_in_time_range(self,timeRangeInMicroSec):
        '''Returns an array of indexes for sample within timeRangeInMicroSec
        timeRangeInMicroSec should be in microseconds'''
        samplesPerRecord = len(self.samples)/self.nRecords
        timeRangeInMicroSec = np.array(np.round(timeRangeInMicroSec),dtype='uint64')
        recordRange = np.searchsorted(self.timestamps,timeRangeInMicroSec)-1
        timeFromRecordStart = (timeRangeInMicroSec-self.timestamps[recordRange])

        sampleInRecord = np.round(1e-6*timeFromRecordStart*self.samplingRate)
        sampleRange = recordRange*samplesPerRecord + sampleInRecord.astype(long)
        sampleIndexes = np.arange(sampleRange[0],sampleRange[-1]+1)
        return sampleIndexes
    def lock_to_event(self,eventOnsetTimes,timeRange):
        '''Make matrix of LFP traces locked to stimulus'''
        if np.any(np.diff(np.diff(self.timestamps))):
            print('Not all LFP records are contiguous. lock_to_event() may not work properly.')
        timeVec = np.arange(timeRange[0],timeRange[-1],1/self.samplingRate)
        nSamples = len(timeVec)
        nTrials = len(eventOnsetTimes)
        lockedLFP = np.empty((nTrials,nSamples))
        for inde,eventTime in enumerate(eventOnsetTimes):
            if not np.isnan(eventTime):
                thisTrialTimeRange = eventTime+timeRange
                sampleIndexes = self.samples_in_time_range(1e6*thisTrialTimeRange)
                ##print '%d %d'%(nSamples,len(sampleIndexes)) ######################### DEBUG
                ####### FIXME: HORRIBLE HACK to avoid errors about sized
                if len(sampleIndexes)!=nSamples:
                    sampleIndexes = sampleIndexes[:-1]
                lockedLFP[inde,:] = self.samples[sampleIndexes]
            else:
                lockedLFP[inde,:] = np.NaN
        return (lockedLFP,timeVec)
    def eventlockedOLD(self,eventOnsetTimes,timeRange):
        print('This function is now OBSOLETE!!! see lock_to_event()')
        wintoplot=np.arange(-20000,40000)
        timeVec = np.arange(timeRange[0],timeRange[-1],1/self.samplingRate)
        firstSample = int(round(timeRange[0]*self.samplingRate))
        nSamples = len(timeVec)
        nTrials = len(eventOnsetTimes)
        samplesRange = np.arange(nSamples,dtype='int')+firstSample
        lockedLFP = np.empty((nTrials,nSamples))
        for inde,eventTime in enumerate(eventOnsetTimes):
            if not np.isnan(eventTime):
                eventTimeInMicroSec = (1e6*eventTime).astype('uint64')
                eventSample = np.searchsorted(self.time,eventTimeInMicroSec)
                thiswin = samplesRange+eventSample
                lockedLFP[inde,:] = self.samples[thiswin]
            else:
                lockedLFP[inde,:] = np.NaN
        return (lockedLFP,timeVec)
    '''
    '''

class DataTetrode(cLoadNeuralynx.DataTetrode):
    '''Access to Neuralynx NTT files containing tetrode data.
    Reading the data is implemented by the parent class from cLoadNeuralynx.'''
    def __init__(self,fileName,readWaves=False):
        cLoadNeuralynx.DataTetrode.__init__(self,fileName,readWaves)
        self.params = {}
        self. _parse_header()
        if readWaves:
            N_CHANNELS = self.params['NumADChannels']
            SAMPLES_PER_SPIKE = self.params['WaveformLength']
            self.samples = self.samples.reshape((N_CHANNELS,SAMPLES_PER_SPIKE,-1),order='F')
    def set_clusters(self,clusterFileOrArray):
        '''Access to KlustaKwik CLU files containing cluster data.'''
        if isinstance(clusterFileOrArray,str):
            self.clusters = np.fromfile(clusterFileOrArray,dtype='int32',sep=' ')[1:]
        else:
            self.clusters = np.array(clusterFileOrArray)
    def _parse_header(self):
        patt = re.compile(r'-NumADChannels (\d+)')
        self.params['NumADChannels'] = int(patt.search(self.header).groups()[0])
        patt = re.compile(r'-WaveformLength (\d+)')
        self.params['WaveformLength'] = int(patt.search(self.header).groups()[0])
        patt = re.compile(r'-SamplingFrequency (\d+\.\d+)')
        self.params['SamplingFrequency'] = float(patt.search(self.header).groups()[0])

class DataEvents(object):
    '''Access to Neuralynx NEV files containing events data
       (pure python implementation)'''
    def __init__(self,fileName,readStrings=False):
        self.filename = fileName
        self._readStrings = readStrings
        self.timestamps = None
        self.eventID = None
        self.valueTTL = None
        self.eventString = None
        self.read_nev(fileName)
        #self.nEvents = self.timestamps.size # Update if read_nev()

    def read_nev(self,fileName):
        HEADERSIZE = 16384 # 16kB
        INFOSIZE = (16+16+16+64 + 5*16 + 32*8)//8
        STRINGSIZE = 128
        RECORDSIZE = INFOSIZE + STRINGSIZE

        fid = open(fileName,'rb')
        header = fid.read(HEADERSIZE)
        # -- Count number of records --
        currentPos = fid.tell()
        fid.seek(0,2) # Go to the end of the file
        fileSize = fid.tell()
        fid.seek(currentPos)
        nRecords = (fileSize-HEADERSIZE)//RECORDSIZE
        self.timestamps = np.empty(nRecords,dtype='uint64')
        self.eventID = np.empty(nRecords,dtype='int16')
        self.valueTTL = np.empty(nRecords,dtype='uint16')
        if self._readStrings:
            self.eventString = np.empty((nRecords,STRINGSIZE),dtype='int8')

        infoFormat = '<hhhQhHhhh' + 8*'l'
        for indrec in range(nRecords):
            oneRecord = fid.read(RECORDSIZE)
            unpackedRecord = unpack(infoFormat,oneRecord[0:INFOSIZE])
            self.timestamps[indrec] = unpackedRecord[3]
            self.eventID[indrec] = unpackedRecord[4]
            self.valueTTL[indrec] = unpackedRecord[5]
            if self._readStrings:
                self.eventString[indrec,:] = np.fromstring(oneRecord[INFOSIZE:],dtype='int8')
        fid.close()

    def test_integrity(self):
        '''Use parity bits to test that there are events for each trial.'''
        # For the following definitions see ExperPort/Settings/Settings_Custom.conf
        bitTRIALIND = 3  # bitID starting from 0
        bitSTIM     = 4
        bitTARGET   = 5 
        bitPARITY1  = 6
        bitPARITY2  = 7
        trialEvents = ( (self.valueTTL & (1<<bitTRIALIND)) != 0 )
        # -- Verify first trials --
        # plot(events.valueTTL[trialEvents]>>6,'.-')
        '''
        clf()
        hold(True)
        plot(trialStartTime  - trialStartTime[1]  ,'.')
        plot(trialStartTimeNL- trialStartTimeNL[0],'o',mfc='none',mec='r')
        hold(False)
        '''
    def find_bit_changes(self,bitInd):
        thisBitArray = (self.valueTTL & (1<<bitInd)) != 0
        valueDiff = np.r_[0,np.diff(thisBitArray.astype(int))]
        onsetInds = valueDiff==1
        offsetInds = valueDiff==-1
        return(onsetInds,offsetInds)
    def __str__(self):
        objStrings = []
        for oneTTL in self.valueTTL:
            objStrings.append('%s\n'%dec2bin(oneTTL,True))
        return ''.join(objStrings)
    def as_string(self,indexes=[],space=True):
        '''Return string representation of TTL values in binary'''
        if not len(indexes):
             indexes = range(self.valueTTL.size)
        for oneTTL in self.valueTTL[indexes]:
            print dec2bin(oneTTL,space)


def dec2bin(n,space=False):
    "Decimal to binary as string"
    n = n+(1<<16)
    bStr = ''
    bitInd = 0
    while n > 0:
        if space and not ((bitInd)%4) and bitInd>0:
            bStr = ' ' + bStr
        bStr = str(n % 2) + bStr
        n = n >> 1
        bitInd += 1
    return bStr


def read_header(dataFile):
    '''Read (only) the header of a neuralynx file.'''
    HEADERSIZE = 16384 # in bytes
    f = open(dataFile, 'r')
    header = f.read(HEADERSIZE)
    f.close()
    return header

'''
def read_clu():
    cfile='/var/data/neuralynx/saja064/2011-02-04_11-28-09_kk/TT2.clu.1'
    fid = open(fileName,'rb')
    time d = np.loadtxt(cfile,dtype='int8')
    time d = np.fromfile(cfile,dtype='int8')
'''

if __name__=='__main__':
    CASE = 3
    if CASE==1:
        data=DataTetrode('/var/tmp/TT0.ntt')
    elif CASE==2:
        '''Load and plot LFP'''
        import pylab as plt
        #dataLFP=DataCont('/var/data/neuralynx/saja125/2012-02-02_17-16-59/CSC17.ncs') # small
        dataLFP=DataCont('/var/data/neuralynx/saja125/2012-02-02_16-33-29/CSC17.ncs') # large
        # lock_to_event
        plt.clf()
        plt.plot(dataLFP.samples[:10000])
        plt.draw()
        plt.show()
    elif CASE==3:
        import pylab as plt
        dataLFP=DataCont('/var/data/neuralynx/saja125/2012-02-02_17-16-59/CSC27.ncs')
        trialOnsetTime = np.array([ 3356.156999,  3361.407959,  3365.967437,  3370.861073])
        #timeRange = [2949053060-33, 2949053060+34]  #2949053060
        #print dataLFP.samples_in_time_range(timeRange)
        (lockedLFP,timeVec) = dataLFP.lock_to_event(trialOnsetTime,[0.38,0.50])
        plt.clf()
        plt.plot(timeVec,lockedLFP.T)
        plt.draw()
        plt.show()

        
