# -*- coding: utf-8 -*-
"""
Receptor interactions in mice

Model the dynamic system (using the control human and mice data), and analyze 
if we can make predictions of interventions, validating with the knockout data.
Nomenclature is from the 2001 Paxinos atlas of the mouse brain


@author: akhan

Step 1: Adjacency matrix A
Step 2: Predict intervention time, predict final intervention state
Step 3: Compare predictions with gene expression

Data:

- Knockout mice_data/densities_5HT1A_KO_model_control_KO_animals.xlsx

"""

import os
import numpy as np
from openpyxl import load_workbook
from scipy.interpolate import griddata

N_META_COLS = 3 # Metadata in Excel file (group, animal, region)

N_MICE_KO = 6
N_MICE_CTRL = 7
N_RECEPTORS = 18
N_REGIONS = 20

path_mice_densities = os.path.join(os.getcwd(), 'data/Knockout mice_data/densities_5HT1A_KO_model_control_KO_animals.xlsx')


# Copy spreadsheet data
wb = load_workbook(filename=path_mice_densities)
ws = wb[wb.sheetnames[0]]

# [group, animal, region, rec1, rec2, ..., recN]

rec_densities = []

# Copy data from spreadsheet
for row in ws.rows:
    r = []
    for cell in row:
        r += [cell.value]
    rec_densities += [r]
    

# Prune empty cells
rec_densities = np.asarray(rec_densities)[: (N_MICE_KO + N_MICE_CTRL)*N_REGIONS + 1, : N_RECEPTORS + N_META_COLS]


# Remove RACL
rec_data = np.delete(rec_densities, (N_META_COLS+N_RECEPTORS-2), axis=1)
N_RECEPTORS -= 1

# No metadata
rec_ctrl = rec_data[1:1+(N_REGIONS*N_MICE_CTRL), N_META_COLS:]
rec_ko = rec_data[1+(N_REGIONS*N_MICE_CTRL):, N_META_COLS:]

# Add another dimension for animals
rec_ctrl = rec_ctrl.reshape((N_MICE_CTRL, N_REGIONS, N_RECEPTORS))
rec_ko = rec_ko.reshape((N_MICE_KO, N_REGIONS, N_RECEPTORS))

# TODO: build sparse adjacency matrices