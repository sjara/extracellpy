#!/usr/env python
'''
This script will create spike feature files, run KlustaKwik and generate reports for one tetrode.
It is useful for running the clustering process in a remote server.

EXAMPLE:
python runclustering.py saja000 2011-04-04_11-54-29 8

'''

import spikesorting
reload(spikesorting)
import sys

__author__ = 'Santiago Jaramillo'
__version__ = '0.1'

#if len(sys.argv)<3: pass  # FIXME: Check number of inputs

animalName   = sys.argv[1]
ephysSession = sys.argv[2]
tetrode      = int(sys.argv[3])

oneTT = spikesorting.TetrodeToCluster(animalName,ephysSession,tetrode)
oneTT.create_fet_files()
oneTT.run_clustering()
oneTT.save_report()
