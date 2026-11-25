# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 18:29:41 2025

@author: meg16
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse
import setisynch.gaia as ss


# Parse the VOTable
votable = parse(r"C:\Users\meg16\Downloads\Data (6)\Data\GCNS_allcol.vot")  # Update with your actual file path
table = votable.get_first_table().to_table()
# Convert to Pandas DataFrame
gaia = table.to_pandas()
print(gaia.columns)

line160 = gaia[-160:]

#%%
planet_data = pd.read_csv(r"C:\Users\meg16\Downloads\Data\Data\habitableworldscatalogue.csv")

spect_filter =  planet_data[planet_data['S_TYPE'].astype(str).str.startswith(('F', 'K', 'G'))]
#%%

line160['Dist50'] = line160['Dist50']*3.6e3
matches = [
    (a, b) for a in spect_filter for b in line160 if a['S_DISTANCE'] == b['Dist50']
]
#%%

# Ensure 'S_DISTANCE' and 'Dist50' are numeric
spect_filter['S_DISTANCE'] = pd.to_numeric(spect_filter['S_DISTANCE'], errors='coerce')
line160['Dist50'] = pd.to_numeric(line160['Dist50'], errors='coerce')
gaia['Dist50'] = pd.to_numeric(gaia['Dist50'], errors='coerce')
# Convert 'Dist50' to the same units (assuming you want light-years)
line160['Dist50'] = line160['Dist50'] * 3.6e3
gaia['Dist50'] = gaia['Dist50'] * 3.6e3

# Find matches based on 'S_DISTANCE' and 'Dist50'
matches = [
    (a, b) for a in spect_filter.to_dict(orient='records')  
    for b in gaia.to_dict(orient='records')  
    if a['S_DISTANCE'] == b['Dist50']
]

