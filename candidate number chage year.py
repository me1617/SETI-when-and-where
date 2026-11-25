# -*- coding: utf-8 -*-
"""
Created on Thu May 15 22:35:08 2025

@author: meg16
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from datetime import datetime, time
import Ellipsoid_Functions as ef
import pickle
from astropy.io.votable import parse
from astropy.coordinates import SkyCoord
from astropy import units as u


print(os.getcwd())
#sys.path.append("/Users/meg16/Downloads/Data (6)/Data/Ellipsoid_Functions.py")



# Parse the VOTable
votable = parse(r"C:\Users\meg16\Downloads\Data (6)\Data\GCNS_allcol.vot")  # Update with your actual file path
table = votable.get_first_table().to_table()

# Convert to Pandas DataFrame
gaia = table.to_pandas()
print(gaia.columns)
#print("Number of stars in GCNS:", len(gaia))
# 331312
# 195978
# 27904
# 767
# Prepare empty lists to hold the indices of rows that meet each condition.
#%%
# Figure 1

SN_dist = 51.4 #kpc
event_name = '1987a'
SN_dist_error = 1.2 #kpc
SN_Dec = -69.2699417
SN_RA = 5.5911111
SN = SkyCoord(SN_RA *u.deg, SN_Dec *u.deg, distance=SN_dist*u.kpc, frame='icrs') 

kpc_to_lyr = 3261.56
c_m_s = 3e8
c_lyr_s = c_m_s/9.461e+15
SN_date = datetime.strptime('1987-02-23 05:31:00', '%Y-%m-%d %H:%M:%S')
# Current date and time
on_ellipsoid_counts = []
for year in range(1988, 2030 + 1):
    
    dt2 = datetime.strptime(f'{year}-02-23 05:31:00', '%Y-%m-%d %H:%M:%S')
    time_since_SN=ef.date_diff_in_Seconds(dt2, SN_date)
    [x_e,y_e]=ef.calc_SETI_ellipsoid(SN_dist*kpc_to_lyr,time_since_SN,'1987a') #lyr
    L = SN_dist*kpc_to_lyr / 2  # lyr
    T = time_since_SN * c_lyr_s  # lyr
    A = (2 * L + T) / 2  # lyr
    
    x_e_centered = x_e-A
    E_pos=[0,0]
    SN_pos=[-51.4,0]
    plt.show()
    seen_SN_indices = []
    in_el_indices = []
    on_el_indices = []
    ellipse_threshold = 0.1
    for i, row in gaia.iterrows():
        # For a given star
        star = SkyCoord(row['RA_ICRS']*u.deg, row['DE_ICRS']*u.deg, row['Dist50']*u.kpc, frame='icrs')
        sep = SN.separation_3d(star)
            
        if ((sep + row['Dist50']*u.kpc).to(u.lyr) >= 2*A*u.lyr - ellipse_threshold*u.lyr) and \
           ((sep + row['Dist50']*u.kpc).to(u.lyr) <= 2*A*u.lyr + ellipse_threshold*u.lyr):
            on_el_indices.append(i)
    on_ellipsoid_counts.append({'Year': year, 'StarsOnEllipse': len(on_el_indices)})
    
    # Now build the new DataFrames using .loc with the accumulated indices.
#     seen_SN = gaia.loc[seen_SN_indices]
#     in_el = gaia.loc[in_el_indices]
#     on_el = gaia.loc[on_el_indices]

# print("Number of stars in GCNS:", len(gaia))
# print("Number of stars seen event:", len(seen_SN))
# print("Number of stars in ellipse:", len(in_el))
# print("Number of stars on ellipse:", len(on_el))


