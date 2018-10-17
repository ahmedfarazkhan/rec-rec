# -*- coding: utf-8 -*-
"""
Map data across parcellations

@author: akhan


Data 
- AtlasesMRI/brodmann.hdr
- AtlasesMRI/brodmann.img
- AtlasesMRI/Julich_Anatomy_v22c.nii
- AtlasesMRI/Julich_Anatomy_v22c_MPM.mat

- Human/receptor_data.txt
- Human/receptor_data.xls
- Human/receptor_data_area_atlas_correlation.xlsx

"""

import os
import numpy as np
import h5py
from openpyxl import load_workbook


N_SUBJECTS = 20 #200
N_REGIONS = 40
N_FACTORS = 6
N_RECEPTORS = 19

def get_MAT_data(field_name, index):
    reference = f_julich_metadata["MAP"][field_name][index][0]
    return f_julich_metadata[reference] 
    
# Match areas  
path_julich_metadata = os.path.join(os.getcwd(), 'data/AtlasesMRI/Julich_Anatomy_v22c_MPM.mat')
path_human_metadata = os.path.join(os.getcwd(), 'data/Human/receptor_data_area_atlas_correlation.xlsx')

###############################################################################

f_julich_metadata = h5py.File(path_julich_metadata, 'r')
julich_map = f_julich_metadata['MAP']

julich_names = f_julich_metadata['MAP']['name']

# List all fields
print("Fields: %s" %f_julich_metadata["MAP"].keys())

field_name = 'name'
julich_area_names = []
for i in range(f_julich_metadata['MAP'][field_name].shape[0]):
    area_name_arr = np.asarray(get_MAT_data(field_name=field_name, index=i)[...])
    julich_area_names += [''.join([chr(c_arr[0]) for c_arr in area_name_arr])]

print julich_area_names

###############################################################################
print "\n"

# Copy spreadsheet data
wb = load_workbook(filename=path_human_metadata)
ws = wb[wb.sheetnames[0]]

# [Region, Julich brain atlas]
area_atlas = []

# Copy data from spreadsheet
for row in ws.rows:#[:1]: 
    r = []
    for cell in row:
        r += [cell.value]
    area_atlas += [r]
    
# Last col is empty
area_atlas = np.delete(area_atlas, (2), axis=1)