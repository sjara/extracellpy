#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Objects and methods for keeping information about isolated cells'''

import numpy as np
import sessionanalysis
reload(sessionanalysis)

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

class CellInfo(object):
    '''Container of information for a given cell. 
    '''
    def __init__(self, animalName='', ephysSession='',behavSession='', tetrode=-1, cluster=-1,
                 trialsToExclude=[]):
        # -- Basic info --
        self.animalName = animalName
        self.ephysSession = ephysSession
        self.behavSession = behavSession
        #self.sessionID = sessionID
        self.tetrode = tetrode
        self.cluster = cluster
        # -- Trial selection --
        self.trialsToExclude = np.array(trialsToExclude,dtype=int)
        # -- Response properties --
        self.soundResponsive = None
    def __repr__(self):
        objStrings = []
        for key,value in sorted(vars(self).iteritems()):
            objStrings.append('%s: %s\n'%(key,str(value)))
        return ''.join(objStrings)
    def __str__(self):
        objStr = '%s %s T%dc%d'%(self.animalName,self.ephysSession,
                                 self.tetrode,self.cluster)
        return objStr
    '''
    def save_locked_spikes(self,outputDir,timeRange=np.array([-0.3,0.9]),lockTo=1):
        cellStr = str(self).replace(' ','_')
        spikesanalysis.save_locked_spikes(self,)
    def info(self):
        # FIXME: Why do I need this? shouldn't __str__ be enough?
        infoStr = '%s [%s] T%dc%d'%(self.animalName,self.behavSession,
                                    self.tetrode,self.cluster)
        return infoStr
    '''

class MultiUnitInfo(object):
    '''Container of information for a multiunit site'''
    def __init__(self, animalName='', ephysSession='',behavSession='', tetrode=-1, clusters=-1):
        '''clusters can be empty (all spikes will be included)'''
        # -- Basic info --
        self.animalName = animalName
        self.ephysSession = ephysSession
        self.behavSession = behavSession
        self.tetrode = tetrode
        self.clusters = clusters
        # -- Response properties --
        #self.soundResponsive = None
    def __repr__(self):
        objStrings = []
        for key,value in sorted(vars(self).iteritems()):
            objStrings.append('%s: %s\n'%(key,str(value)))
        return ''.join(objStrings)
    def __str__(self):
        objStr = '%s %s T%d'%(self.animalName,self.ephysSession,
                              self.tetrode)
        return objStr
    
class LFPInfo(object):
    '''Container of information for a given LFP. 
    '''
    def __init__(self, animalName='', ephysSession='',behavSession='', electrode=-1):
        # -- Basic info --
        self.animalName = animalName
        self.ephysSession = ephysSession
        self.behavSession = behavSession
        self.electrode = electrode
    def __repr__(self):
        objStrings = []
        for key,value in sorted(vars(self).iteritems()):
            objStrings.append('%s: %s\n'%(key,str(value)))
        return ''.join(objStrings)
    def __str__(self):
        objStr = '%s %s E%d'%(self.animalName,self.ephysSession,
                                 self.electrode)
        return objStr


class CellDatabase(list):
    '''Container of set of cells.
    '''
    def __init__(self):
        super(CellDatabase, self).__init__()
    #def extend(self,extraCellDB):  ### NOT NECESSARY
    #    list.extend(self,extraCellDB)
    def findcell(self,firstParam,behavSession='',tetrode=-1,cluster=-1):
        '''Find index of cell. It can be used in two ways:
        cellDB.findcell('saja099','20110528a',1,11)
        cellDB.findcell(onecell)
        '''
        if isinstance(firstParam,str):
            onecell = CellInfo(firstParam,'',behavSession,tetrode,cluster)
        else:
            onecell = firstParam
        cellIndex = None
        for ind,cell in enumerate(self):
            if onecell.animalName==cell.animalName:
                if onecell.behavSession==cell.behavSession:
                    if onecell.tetrode==cell.tetrode:
                        if onecell.cluster==cell.cluster:
                            cellIndex = ind
        return cellIndex
    def set_soundResponsive(self,zScores,threshold=3):
        '''Set soundResponsive flag for each cell, given zScores
        
        zScores: numpy array (nTimeBins,nConditions,nCells)
        threshold: above this or below negative this it is considered responsive
        '''
        for indcell,onecell in enumerate(self):
            onecell.soundResponsive = np.any(abs(zScores[:,:,indcell])>threshold)
    def get_vector(self,varname):
        '''EXAMPLE: cellDB.get_vector('tetrode')
        '''
        return np.array([getattr(onecell, varname) for onecell in self])
    def subset(self,indexes):
        subsetDB = CellDatabase()
        for ind in indexes:
            subsetDB.append(self[ind])
        return subsetDB
    def __str__(self):
        objStrings = []
        for ind,c in enumerate(self):
            objStrings.append('[%d] %s\n'%(ind,c))
        return ''.join(objStrings)
    def save_locked_spikes(self,outputDir,timeRange=np.array([-0.3,0.9]),lockTo=1):
        sessionanalysis.save_data_each_cell(self,outputDir,timeRange=timeRange,lockTo=lockTo)
    def evaluate_response(self):
        from extracellpy import sessionanalysis
        ############# FINISH THIS #########
    def evaluate_responsiveness(self):
        pass


class MultiUnitDatabase(list):
    '''Container of set of multiunit sites.
    '''
    def __init__(self):
        super(MultiUnitDatabase, self).__init__()
    def __str__(self):
        objStrings = []
        for ind,c in enumerate(self):
            objStrings.append('[%d] %s\n'%(ind,c))
        return ''.join(objStrings)
    def save_locked_spikes(self,outputDir,timeRange=np.array([-0.3,0.9]),lockTo=1):
        sessionanalysis.save_data_each_mu(self,outputDir,timeRange=timeRange,lockTo=1)



'''
run -e allcells_saja099.py
alldata=np.load('/var/tmp/saja099_responsiveness.npz')
cellDB.set_soundResponsive(alldata['zStatsEachCell'])
'''

