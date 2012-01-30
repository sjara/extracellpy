#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Functions for loading spikes, clusters and LFP
'''

import numpy
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


def read_t(SpikesDataFileName):
    '''Reads data from an MClust .t file.
       These files contain events/spikes data.
       Based on MClust-3.4/Matlab/SpecificUtils/ReadHeader.m
    '''
    fid = open(SpikesDataFileName,'rb')
    (DataPosition,HeaderText) = read_t_header(fid)
    fid.seek(DataPosition)
    ByteString = fid.read()
    Nsamples = len(ByteString)/4
    timestamps = numpy.empty(Nsamples,dtype=numpy.uint32)
    for inds in range(Nsamples):
        timestamps[inds] = unpack('>I', ByteString[4*inds:4*inds+4])[0]
    #timestamps = numpy.fromfile(fid,dtype='uint32')
    fid.close()
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
    timestamps = numpy.fromfile(fid,dtype='uint32')[0::nBytesPerWaverform/32+1]
    fid.close()
    timestamps = 1e-4*timestamps        # From 100usec to sec (see NSpike/MClust docs)
    return timestamps

