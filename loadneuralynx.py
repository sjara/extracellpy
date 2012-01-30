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
        self.__test_integrity()
        patt = re.compile(r'-SamplingFrequency (\d+\.\d+)')
        self.samplingRate = float(patt.search(self.header).groups()[0])
        self.firstTimestamp = self.timestamps[0]
        self.time = self.generate_timevec()
    def __test_integrity(self):
        if np.any(np.diff(np.diff(self.timestamps))):
            raise TypeError('Not all records are contiguous. Packets lost?')
    def generate_timevec(self):
        '''This may be a more complicated than necessary: I'm creating an equally
        spaced vector between each record timestamp.
        '''
        dRec = self.timestamps[1]-self.timestamps[0]
        samplesPerRecord = len(self.samples)/self.nRecords
        timeArray = np.tile(self.timestamps.reshape(self.timestamps.shape[0], -1),
                            (1,samplesPerRecord))
        interimTimestamps = np.linspace(0,dRec,samplesPerRecord,endpoint=False).astype(timeArray.dtype)
        timeVector = timeArray+interimTimestamps
        return timeVector.flatten()
    def __parse_header(self):
        pass
    #def __repr__(self):
    #    pass
        '''
        objStrings = []
        for key,value in sorted(vars(self).iteritems()):
            objStrings.append('%s: %s\n'%(key,str(value)))
        return ''.join(objStrings)
        '''
    def eventlocked(self,eventOnsetTimes,timeRange):
        '''Make matrix of LFP traces locked to stimulus'''
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

class DataTetrode(cLoadNeuralynx.DataTetrode):
    '''Access to Neuralynx NTT files containing tetrode data.
    Reading the data is implemented by the parent class from cLoadNeuralynx.'''
    def __init__(self,fileName,readWaves=False):
        cLoadNeuralynx.DataTetrode.__init__(self,fileName,readWaves)
    def set_clusters(self,clusterFile):
        '''Access to KlustaKwik CLU files containing cluster data.'''
        self.clusters = np.fromfile(clusterFile,dtype='int32',sep=' ')[1:]


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


'''
def read_clu():
    cfile='/var/data/neuralynx/saja064/2011-02-04_11-28-09_kk/TT2.clu.1'
    fid = open(fileName,'rb')
    time d = np.loadtxt(cfile,dtype='int8')
    time d = np.fromfile(cfile,dtype='int8')
'''

if __name__=='__main__':
    data=DataTetrode('/var/tmp/TT0.ntt')
