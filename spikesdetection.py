#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Detection of spikes from continuous data.
'''

from extracellpy import settings
import numpy as np
import os

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

EPHYSPATH = settings.EPHYS_PATH

def filter_continuous(data,samplingRate,highPassCut=600,order=6):
    '''
    Filter continuous signal to extract spikes
    '''
    import scipy.signal as signal
    Wn = highPassCut*2/samplingRate # Normalized
    [b,a]=signal.butter(order, Wn, btype='high')
    print('Filtering continuous data...')
    filteredData = signal.lfilter(b,a,data)
    print('Done!')
    return filteredData


def find_threshold_crossing(data,threshold=None,refractorySamples=8,dual=False):
    '''dual defines if detection is done for both positive and negative crossings'''
    if not threshold:
        threshold = np.std(data)
    crossesBool = data>threshold
    if dual:
        crossesBool = crossesBool | (data<-threshold)
    crosses = np.flatnonzero(crossesBool)
    beyondRefraction = np.diff(np.r_[0,crosses])>refractorySamples
    crossesCleaned = crosses[beyondRefraction]
    params = {'threshold':threshold,'refractorySamples':refractorySamples,
              'dual':dual}
    return crossesCleaned,params

class ExtractedSpikes(object):
    #def __init__(self):
    def __init__(self,spikesFileName=None,contFileName=None):
        self.timestamps = None
        self.detectionParams = {}
        if spikesFileName:
            try:
                self.from_file(spikesFileName)
            except IOError:
                # TODO: Ask if want to generate
                raise
    def from_file(self,fileName):
        import h5py
        # TODO: read also 'header' info.
        dataFile = h5py.File(fileName,'r')
        self.timestamps = dataFile['timestamps'][...]
        dataFile.close()
    def extract_spikes(self,animalName,ephysSession,electrode,
                       threshold=None,refractorySamples=8,dual=False,showSignal=False):
        import h5py
        import loadneuralynx
        # -- Load continuous data --
        dataDir = os.path.join(EPHYSPATH,'%s/%s/'%(animalName,ephysSession))
        contDataFile = os.path.join(dataDir,'CSC%d.ncs'%electrode)
        dataLFP = loadneuralynx.DataCont(contDataFile)
        # -- Extract spikes --
        filtLFP = filter_continuous(dataLFP.samples,dataLFP.samplingRate)
        params = {'threshold':threshold,'refractorySamples':refractorySamples,
                  'dual':dual}
        (crosses,params) = find_threshold_crossing(filtLFP,**params)
        spikeTimes = dataLFP.time_of_sample(crosses)
        self.detectionParams = params
        self.timestamps = spikeTimes
        # -- Save spikes --
        spikesDataDir = settings.EXTRACTED_SPIKES_PATH%animalName
        spikesFileName = '%s_%s_e%02d_spikes.h5'%(animalName,ephysSession,electrode)
        spikesFileFull = os.path.join(spikesDataDir,spikesFileName)
        print 'Saving data to %s'%spikesFileFull
        dataFile = h5py.File(spikesFileFull,'w')
        try:
            dataFile.create_dataset('timestamps', data=self.timestamps)
        except:
            dataFile.close() # Close before raising error
            raise
        dataFile.close()
        if showSignal:
            import pylab as plt
            plt.clf()
            plt.plot(filtLFP[:10000],'b-')
            plt.hold(True)
            plt.axhline(threshold,color='0.75')
            #plt.plot(crosses[crosses<10000],2*np.tile(threshold,len(crosses)),'.k')
            plt.hold(False)
            plt.draw()
            plt.show()
        '''
        '''
        return spikesFileFull
    
if __name__ == "__main__":
    CASE = 4
    if CASE==1:
        from pylab import *
        import loadneuralynx
        import scipy.signal as signal
        #import lfilter, bessel

        fileName = '/var/data/neuralynx/saja125/2012-02-03_18-32-48/CSC8.ncs' # Large spikes
        #fileName = '/var/data/neuralynx/saja125/2012-02-05_17-34-15/CSC2.ncs' # Very large spikes
        dataLFP = loadneuralynx.DataCont(fileName)

        # -- Filter signal (high-pass) --
        highpassFreq = 100 # Hz
        Wn = highpassFreq*2/dataLFP.samplingRate # Normalized
        [b,a]=signal.butter(6, Wn, btype='high')
        filtLFP = signal.lfilter(b,a,dataLFP.samples[:500000])
        #filtLFP = signal.filtfilt(b,a,dataLFP.samples[:500000])

        highpassFreq = 600 # Hz
        Wn = highpassFreq*2/dataLFP.samplingRate # Normalized
        [b,a]=signal.butter(6, Wn, btype='high')
        filtLFP2 = signal.lfilter(b,a,dataLFP.samples[:500000])


        rangeToPlot = np.arange(100000)+200000
        plot(dataLFP.samples[rangeToPlot],'-')
        hold(True)
        plot(filtLFP[rangeToPlot],'g-')
        plot(filtLFP2[rangeToPlot],'r-')
        hold(False)
        #xlim([72000,84000])
        #xlim([90000,92000])
        #xlim([85400,86600]) # Good for 2012-02-03_18-32-48/CSC8.ncs'
        xlim([82300,83000]) # Good for 2012-02-03_18-32-48/CSC8.ncs'
        ylim([-6000,8000])
        draw()
        show()

    if CASE==2:
        from pylab import *
        import loadneuralynx
        reload(loadneuralynx)
        import scipy.signal as signal

        #import line_profiler as lprof
        #profile = lprof.LineProfiler(find_threshold_crossing)
        #profile.enable()

        fileName = '/var/data/neuralynx/saja125/2012-02-03_18-32-48/CSC8.ncs' # Large spikes
        dataLFP = loadneuralynx.DataCont(fileName)
        dataChunk = dataLFP.samples[arange(10000)+200000]
        #dataChunk = dataLFP.samples

        # -- Filter signal (high-pass) --
        filtLFP = filter_continuous(dataChunk,dataLFP.samplingRate)

        # -- Find crossings --
        threshold = 2000
        refractorySamples = 8
        crosses,params = find_threshold_crossing(filtLFP,threshold,refractorySamples,dual=True)
        
        spikeTimes = dataLFP.time_of_sample(crosses)
        #profile.print_stats()

        #rangeToPlot = np.arange(10000)
        plot(filtLFP,'b-')
        hold(True)
        axhline(threshold,color='0.75')
        plot(crosses,2*tile(threshold,len(crosses)),'.k')
        hold(False)
        draw()
        show()
        '''
        '''
    if CASE==3:
        animalName = 'saja125'
        ephysSession = '2012-02-03_18-32-48'
        electrode = 8
        spikes = ExtractedSpikes()
        dataFile = spikes.extract_spikes(animalName,ephysSession,electrode,threshold=2000,dual=True)
        #newspikes = ExtractedSpikes(dataFile)
    if CASE==4:
        from pylab import *
        animalName = 'saja125'
        ephysSession = '2012-02-03_18-32-48'
        electrode = 8
        spikesDataDir = settings.EXTRACTED_SPIKES_PATH%animalName
        spikesFileName = '%s_%s_e%02d_spikes.h5'%(animalName,ephysSession,electrode)
        spikesFileFull = os.path.join(spikesDataDir,spikesFileName)
        newspikes = ExtractedSpikes(spikesFileFull)
        
        '''
        # -- Filter signal (high-pass) --
        filtLFP = filter_continuous(dataChunk,dataLFP.samplingRate)

        # -- Find crossings --
        threshold = 2000
        refractorySamples = 8
        crosses = find_threshold_crossing(filtLFP,threshold,refractorySamples,dual=True)
        
        spikeTimes = dataLFP.time_of_sample(crosses)
        #profile.print_stats()

        #rangeToPlot = np.arange(10000)
        plot(filtLFP,'b-')
        hold(True)
        axhline(threshold,color='0.75')
        plot(crosses,2*tile(threshold,len(crosses)),'.k')
        hold(False)
        draw()
        show()
        '''

