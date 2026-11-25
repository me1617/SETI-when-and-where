# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 17:43:47 2025

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
import setisynch.gaia as ss
import setisynch.hygdata_v3 as hs



# Open all the files that come from the GCNS data for the threshold value of 0.5
planet_data_seen = pd.read_csv(r"C:\Users\meg16\Downloads\Data (6)\Data\seen_SN.csv")

planet_data_on = pd.read_csv(r"C:\Users\meg16\Downloads\Data (6)\Data\on_el_th0.5lyr.csv")

planet_data_in = pd.read_csv(r"C:\Users\meg16\Downloads\Data (6)\Data\in_el.csv")

votable = parse(r"C:\Users\meg16\Downloads\Data (6)\Data\GCNS_allcol.vot")  
table = votable.get_first_table().to_table()

# Convert to Pandas DataFrame
gaia = table.to_pandas()
print(gaia.columns)
# Filter for stars visible by the Green Bank Telescope
visible_planets_seen = planet_data_seen[(planet_data_seen['Dist50']*3.26156e3<=100000) & (planet_data_seen['DE_ICRS'] >= -46) & (planet_data_seen['DE_ICRS'] <= 90) &
                                        (planet_data_seen['RA_ICRS'] >= 0) & (planet_data_seen['RA_ICRS'] <=24)]

visible_planets_on = planet_data_on[(planet_data_on['Dist50']*3.26156e3<=100000) & (planet_data_on['DE_ICRS'] >= -46) & (planet_data_on['DE_ICRS'] <= 90)  &
                                    (planet_data_on['RA_ICRS'] >= 0) & (planet_data_on['RA_ICRS'] <=24)]

visible_planets_in = planet_data_in[(planet_data_in['Dist50']*3.26156e3<=100000) & (planet_data_in['DE_ICRS'] >= -46) & (planet_data_in['DE_ICRS'] <= 90) &
                                    (planet_data_in['RA_ICRS'] >= 0) & (planet_data_in['RA_ICRS'] <=24)]

visible_planets_GCNS = gaia

# Obtaining the G magnitude of the datasets and the BP and RP magnitudes
gmag_seen = visible_planets_seen['gmag']
bpmag_seen = visible_planets_seen['BPmag']
rpmag_seen = visible_planets_seen['RPmag']

gmag_GCNS = visible_planets_GCNS['gmag']
bpmag_GCNS = visible_planets_GCNS['BPmag']
rpmag_GCNS = visible_planets_GCNS['RPmag']

gmag_on = visible_planets_on['gmag']
bpmag_on = visible_planets_on['BPmag']
rpmag_on = visible_planets_on['RPmag']

gmag_in = visible_planets_in['gmag']
bpmag_in = visible_planets_in['BPmag']
rpmag_in = visible_planets_in['RPmag']

# Calculate the x axis of each file by subtracting the bp against the rp
x_mag_seen = bpmag_seen-rpmag_seen

x_mag_on = bpmag_on-rpmag_on

x_mag_in = bpmag_in-rpmag_in

x_mag_GCNS = bpmag_GCNS-rpmag_GCNS

#%%

import matplotlib.cm as cm
# Combine the data into one array to plot easier for the colourmap
x_mag = np.concatenate([ x_mag_seen, x_mag_in, x_mag_on])
gmag = np.concatenate([ gmag_seen, gmag_in, gmag_on])
color_var = np.concatenate([ gmag_seen, gmag_in, gmag_on]) 

# Create the scatter plot with the color mapped according to the continuous variable
fig= plt.subplots(figsize=(12, 5))
sc = plt.scatter(x_mag, gmag, c=color_var, cmap=cm.plasma, s=10)
plt.xlabel('bp-rp Magnitude')
plt.ylabel('Gmag')
plt.show()

fig= plt.subplots(figsize=(12, 5))
#plt.scatter(x_mag_GCNS, gmag_GCNS, color = 'orange', s=    10)
plt.scatter(x_mag_seen, gmag_seen, color = 'green', s= 10)

plt.scatter(x_mag_on, gmag_on, color = 'red', s= 10)
plt.scatter(x_mag_in, gmag_in, color = 'blue', s= 10)
plt.xlabel('bp-rp Magnitude')
plt.ylabel('Gmag')
plt.show()
#%%

x_mag_GCNS_b = (bpmag_GCNS-rpmag_GCNS)<0
x_mag_GCNS_bw = ((bpmag_GCNS-rpmag_GCNS)>0) & ((bpmag_GCNS-rpmag_GCNS)<0.5)
x_mag_GCNS_wy = ((bpmag_GCNS-rpmag_GCNS)>0.5) & ((bpmag_GCNS-rpmag_GCNS)<1.5)
x_mag_GCNS_or = ((bpmag_GCNS-rpmag_GCNS)>1.5) & ((bpmag_GCNS-rpmag_GCNS)<3)
x_mag_GCNS_r = (bpmag_GCNS-rpmag_GCNS)>3

#plt.scatter([x_mag_GCNS_b,x_mag_GCNS_bw,x_mag_GCNS_wy,x_mag_GCNS_or,x_mag_GCNS_r],  gmag_GCNS, color = 'orange')
plt.figure(figsize=(8, 6))

# Color-magnitude diagram: (BP - RP) on x-axis, Gmag on y-axis (reversed)
plt.scatter((bpmag_GCNS - rpmag_GCNS)[x_mag_GCNS_b],  gmag_GCNS[x_mag_GCNS_b],  color='blue',   label='Very Blue')
plt.scatter((bpmag_GCNS - rpmag_GCNS)[x_mag_GCNS_bw], gmag_GCNS[x_mag_GCNS_bw], color='cyan',   label='Blue-White')
plt.scatter((bpmag_GCNS - rpmag_GCNS)[x_mag_GCNS_wy], gmag_GCNS[x_mag_GCNS_wy], color='yellow', label='White-Yellow')
plt.scatter((bpmag_GCNS - rpmag_GCNS)[x_mag_GCNS_or], gmag_GCNS[x_mag_GCNS_or], color='orange', label='Orange-Red')
plt.scatter((bpmag_GCNS - rpmag_GCNS)[x_mag_GCNS_r],  gmag_GCNS[x_mag_GCNS_r],  color='red',    label='Very Red')

plt.gca().invert_yaxis()  # Bright stars on top
plt.xlabel('BP - RP (Color Index)')
plt.ylabel('Gmag (Absolute Magnitude)')
plt.title('Color-Magnitude Diagram (Gaia GCNS)')
plt.legend()
#plt.grid(True)
plt.show()

#%%

x_mag_seen_b = (bpmag_seen-rpmag_seen)<0
x_mag_seen_bw = ((bpmag_seen-rpmag_seen)>0) & ((bpmag_seen-rpmag_seen)<0.5)
x_mag_seen_wy = ((bpmag_seen-rpmag_seen)>0.5) & ((bpmag_seen-rpmag_seen)<1.5)
x_mag_seen_or = ((bpmag_seen-rpmag_seen)>1.5) & ((bpmag_seen-rpmag_seen)<3)
x_mag_seen_r = (bpmag_seen-rpmag_seen)>3

x_mag_in_b = (bpmag_in-rpmag_in)<0
x_mag_in_bw = ((bpmag_in-rpmag_in)>0) & ((bpmag_in-rpmag_in)<0.5)
x_mag_in_wy = ((bpmag_in-rpmag_in)>0.5) & ((bpmag_in-rpmag_in)<1.5)
x_mag_in_or = ((bpmag_in-rpmag_in)>1.5) & ((bpmag_in-rpmag_in)<3)
x_mag_in_r = (bpmag_in-rpmag_in)>3

x_mag_on_b = (bpmag_on-rpmag_on)<0
x_mag_on_bw = ((bpmag_on-rpmag_on)>0) & ((bpmag_on-rpmag_on)<0.5)
x_mag_on_wy = ((bpmag_on-rpmag_on)>0.5) & ((bpmag_on-rpmag_on)<1.5)
x_mag_on_or = ((bpmag_on-rpmag_on)>1.5) & ((bpmag_on-rpmag_on)<3)
x_mag_on_r = (bpmag_on-rpmag_on)>3

plt.scatter(x_mag_seen_b,  gmag_seen,  color='blue',   label='Very Blue')
plt.scatter(x_mag_seen_bw,  gmag_seen, color='cyan',   label='Blue-White')
plt.scatter(x_mag_seen_wy,  gmag_seen, color='yellow', label='White-Yellow')
plt.scatter(x_mag_seen_or,  gmag_seen, color='orange', label='Orange-Red')
plt.scatter(x_mag_seen_r,  gmag_seen,   color='red',    label='Very Red')

plt.scatter(x_mag_in_b,  gmag_in,  color='blue',   label='Very Blue')
plt.scatter(x_mag_in_bw,  gmag_in, color='cyan',   label='Blue-White')
plt.scatter(x_mag_in_wy,  gmag_in, color='yellow', label='White-Yellow')
plt.scatter(x_mag_in_or,  gmag_in, color='orange', label='Orange-Red')
plt.scatter(x_mag_in_r,  gmag_in,   color='red',    label='Very Red')

plt.scatter(x_mag_on_b,  gmag_on,  color='blue',   label='Very Blue')
plt.scatter(x_mag_on_bw,  gmag_on, color='cyan',   label='Blue-White')
plt.scatter(x_mag_on_wy,  gmag_on, color='yellow', label='White-Yellow')
plt.scatter(x_mag_on_or,  gmag_on, color='orange', label='Orange-Red')
plt.scatter(x_mag_on_r,  gmag_on,   color='red',    label='Very Red')




