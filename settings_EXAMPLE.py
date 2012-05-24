# -*- coding: utf-8 -*-
'''
General settings for extracellpy.

The string %s will be replaced by the animal's name, for example:
'/tmp/%s_processed/' becomes '/tmp/saja000_processed/'

'''

import os

BEHAVIOR_PATH = '/var/data/BControlData/Data/santiago'
EPHYS_PATH = '/var/data/neuralynx'

###PROCESSED_PATH = os.path.join(EPHYS_PATH)

EXTRACTED_SPIKES_PATH = '/var/data/neuralynx/%s_processed/extractedspikes'

# --- WORKFLOW REVERSAL ---
# NOTE: the symbol %s needs to be there. The scripts will replace it
#       by the animal name
CELL_LIST_PATH = '/home/sjara/zadorlab/data_analysis/extracellpytest'
PROCESSED_REVERSAL_PATH = '/var/data/neuralynx/%s_test/'
#PROCESSED_REVERSAL_PATH = '/var/data/neuralynx/%s_processed/'
#PROCESSED_REVERSAL_PATH = '/var/data/neuralynx/%s_temp/'

###PROCESSED_REVERSAL_PATH = '/var/data/neuralynx/%s_processed/'
###RASTER_PLOTS_REVERSAL_PATH = '/var/data/neuralynx/%s_rasters/'

# -- DO NOT EDIT BELOW ---
CLUSTERS_REPORTS_DIR = 'cluster_reports'

