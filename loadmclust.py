#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Functions for loading spikes, clusters and LFP
'''

import numpy as np
from struct import unpack

__author__ = 'Santiago Jaramillo'
__version__ = "$Revision: 2009-03-06$"


def read_t_header(fid):
    '''Read the header of an MClust .t file
        Example: /var/data/nspike/saja034/saja034_20081205-01spike_mclust/05-000_1.t
       TO DO: check beginheader and report if there is no header (some files may not have it)
    '''
    beginheader = '%%BEGINHEADER\n';
    endheader = '%%ENDHEADER\n';
    HeaderText = ''
    while True:
        line=fid.readline()
        HeaderText += line
        if line==endheader:
            break
    DataPosition = fid.tell()
    return (DataPosition,HeaderText)


def read_t(SpikesDataFileName,insec=True):
    '''Reads data from an MClust .t file.
       These files contain events/spikes data.
       Based on MClust-3.4/Matlab/SpecificUtils/ReadHeader.m

       If 'insec' is True, timestamps are given in seconds (float)
       otherwise they are given as tenths of msec (int32)
    '''
    fid = open(SpikesDataFileName,'rb')
    (DataPosition,HeaderText) = read_t_header(fid)
    fid.seek(DataPosition)
    ByteString = fid.read()
    Nsamples = len(ByteString)/4
    timestamps = np.empty(Nsamples,dtype=np.uint32)
    for inds in range(Nsamples):
        timestamps[inds] = unpack('>I', ByteString[4*inds:4*inds+4])[0]
    #timestamps = np.fromfile(fid,dtype='uint32')
    fid.close()
    if insec:
        timestamps = 1e-4*timestamps        # From 100usec to sec (see NSpike/MClust docs)
    return timestamps

def read_tt(SpikesDataFileName):
    '''Reads data from an NSpike .tt file.
       These files contain events/spikes data.

    FINISH THIS DOC (see test02_loadspikes.py)

% Spike data consists of records of the form: 
%     timestamp(uint32) 
%     waveform(4 chans * 40 points = 160 int16)
    '''
    fid = open(SpikesDataFileName,'rb')
    junk = fid.read(200)
    headersize = junk.find('ENDHEADER')+10
    fid.seek(headersize,0)
    nBytesPerWaverform = 4*40*16
    timestamps = np.fromfile(fid,dtype='uint32')[0::nBytesPerWaverform/32+1]
    fid.close()
    timestamps = 1e-4*timestamps        # From 100usec to sec (see NSpike/MClust docs)
    return timestamps

# ============== From MClust t-files to KK clu-files ==================
'''
We need a way to get data saved by MClust and load it into our processing workflow.
It means the data should be located in a standard folder, named in a particular way
and can be referenced from an allcells python file as:
    if SESSIONID==9999:
        animalName   = 'saja100'
        ephysSession = '2011-11-27_14-49-59'
        behavSession = '20111127a'
        clustersEachTetrode = {7:[8]}
Corresponding to /var/data/neuralynx/saja100/2011-11-27_14-49-59_kk/TT7.clu.1

What we can do is to create a new clu files, for example TT107.clu.1 that contains
the spikes we got from MClust.

The output from MClust is a tfile which contains:
  the timestamps of spikes of interest
The clu file contains:
  the index of the cluster each spike belongs to


'''
def tfiles_to_clu(tFiles,tetrodeFile,outputFilename=''):
    '''Create a clu file from MClust t-files.'''
    from extracellpy import loadneuralynx
    if isinstance(tFiles,str):
        tFiles=[tFiles]
    dataTT = loadneuralynx.DataTetrode(tetrodeFile)
    ###nlts = (dataTT.timestamps/100).astype('int32')
    # Necessary because rounding is different in numpy and matlab
    #  searchsorted (below) will find the right index.
    nlts = np.ceil(dataTT.timestamps.astype(float)/100).astype('int32')
    clusterAssignment = np.zeros(nlts.shape,dtype='int')
    for indfile,tfilename in enumerate(tFiles):
        print 'Processing t-file #%d ...'%(indfile+1)
        tempClusterAssig = np.zeros(nlts.shape,dtype='int')
        ts = read_t(tfilename,insec=False)
        sorted_index = np.searchsorted(nlts,ts,'left')
        clusterAssignment[sorted_index]=(indfile+1)
        '''
        for ind,onets in enumerate(ts):
            spikeMatch = abs(nlts-onets)<2 # Needed because of rounding issue
            # NOTE: this is faster than using flatnonzero and assigning that index
            tempClusterAssig = tempClusterAssig + spikeMatch
        clusterAssignment = clusterAssignment+(indfile+1)*tempClusterAssig
        '''
        ###print ts[:4]
        ###print dataTT.timestamps[clusterAssignment==(indfile+1)][:4]
    if outputFilename:
        print 'Saving clu file to %s'%outputFilename
        dataToSave = np.concatenate(([len(tFiles)+1],clusterAssignment))
        np.savetxt(outputFilename, dataToSave, fmt="%d")
    # -- Save data --
    return clusterAssignment


if __name__ == "__main__":
    CASE = 2
    if CASE==1:
        from extracellpy import loadneuralynx
        tfilename = '/home/sjara/tmp/mclustdata/TT7_2.t'
        ts = read_t(tfilename,insec=False)
        #ts = read_t(tfilename,insec=True)
        tetrodeFile = '/home/sjara/tmp/mclustdata/TT7.ntt'
        dataTT = loadneuralynx.DataTetrode(tetrodeFile)
        ###dataTT.timestamps = dataTT.timestamps.astype(np.float64)*1e-6  # in sec
        #nlts = (dataTT.timestamps/100).astype('int32')
        #nlts = np.round(dataTT.timestamps.astype(float)/100).astype('int32')
        nlts = np.ceil(dataTT.timestamps.astype(float)/100).astype('int32')
        # NOTE: nlts[45] a many others have to rounded to match MClust value
        # NOTE: compare ts[163] and nlts[1202]. MClust and numpy round numbers differently!
        #       dataTT.timestamps[1202]=409680450    ts[163]=4096805
        
        clusterAssignment = np.zeros(nlts.shape,dtype='int')
        sorted_index = np.searchsorted(nlts,ts,'left')
        clusterAssignment[sorted_index]=1
        '''
        TEST:
        nlts[clusterAssignment.astype(bool)]
        ts
        '''
        
        '''
        for ind,onets in enumerate(ts):
            #spikeMatch = abs(nlts-onets)<2 # Needed because of rounding issue
            #valdiff = nlts-onets
            spikeMatch = nlts==onets
            #spikeMatch = (nlts==onets)
            #if not any(spikeMatch): break
            clusterAssignment = clusterAssignment + spikeMatch
            # NOTE: faster than using flatnonzero and assigning that index
        '''
        '''
        Verify sizes:
        len(ts)
        sum(clusterAssignment)
        '''
    elif CASE==2:
        tFiles = ['/home/sjara/tmp/mclustdata/TT7_1.t','/home/sjara/tmp/mclustdata/TT7_2.t']
        tetrodeFile = '/home/sjara/tmp/mclustdata/TT7.ntt'
        outputFilename = '/tmp/TT7.clu.1'
        clu = tfiles_to_clu(tFiles,tetrodeFile,outputFilename)
    elif CASE==3:
        tFiles = ['/home/sjara/tmp/mclustdata/TT7_2.t']
        tetrodeFile = '/home/sjara/tmp/mclustdata/TT7.ntt'
        clu = tfiles_to_clu(tFiles,tetrodeFile)

'''
import numpy as np
x = np.array([3,5,7,1,9,8,6,6])
y = np.array([2,1,5,10,100,6])

index = np.argsort(x)
sorted_x = x[index]
sorted_index = np.searchsorted(sorted_x, y)

yindex = np.take(index, sorted_index, mode="clip")
mask = x[yindex] != y

result = np.ma.array(yindex, mask=mask)
print result
[-- 3 1 -- -- 6]
'''
