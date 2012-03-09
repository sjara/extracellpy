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
