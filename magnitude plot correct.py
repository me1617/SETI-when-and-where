# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 13:45:04 2025

@author: meg16
"""


import pandas as pd
import matplotlib.pyplot as plt
from astropy.io.votable import parse
import matplotlib.cm as cm



#Opening all the GCNS files, including the in, on and seen
planet_data_seen = pd.read_csv(r"C:\Users\meg16\Downloads\Data (6)\Data\seen_SN.csv")

planet_data_on = pd.read_csv(r"C:\Users\meg16\Downloads\Data (6)\Data\on_el_th0.5lyr.csv")

planet_data_in = pd.read_csv(r"C:\Users\meg16\Downloads\Data (6)\Data\in_el.csv")

votable = parse(r"C:\Users\meg16\Downloads\Data (6)\Data\GCNS_allcol.vot") # opening GCNS as vot 
table = votable.get_first_table().to_table() # converting vot to table

gaia = table.to_pandas() # converting table to panda format 
print(gaia.columns)


# Extracting the distances for all data files and coverting from kpc to light years
planet_dist_seen = planet_data_seen['Dist50']*3.26156e3
planet_dist_on = planet_data_on['Dist50']*3.26156e3
planet_dist_in = planet_data_in['Dist50']*3.26156e3
planet_dist_GCNS = gaia['Dist50']*3.26156e3

# Extracting all the magnitudes from G magnitude to the bp and rp magnitudes for each dataset
gmag_seen = planet_data_seen['gmag']
bpmag_seen = planet_data_seen['BPmag']
rpmag_seen = planet_data_seen['RPmag']

gmag_GCNS = gaia['gmag']
bpmag_GCNS = gaia['BPmag']
rpmag_GCNS = gaia['RPmag']

gmag_on = planet_data_on['gmag']
bpmag_on = planet_data_on['BPmag']
rpmag_on = planet_data_on['RPmag']

gmag_in = planet_data_in['gmag']
bpmag_in = planet_data_in['BPmag']
rpmag_in = planet_data_in['RPmag']

#Calculating the x variable from the bp-rp calculation to create the spectra
x_mag_seen = bpmag_seen-rpmag_seen

x_mag_on = bpmag_on-rpmag_on

x_mag_in = bpmag_in-rpmag_in

x_mag_GCNS = bpmag_GCNS-rpmag_GCNS

# x and y axis for the Sun based off the magnitudes
y_sun = 4.68 # Gmag

x_sun = -0.13+2.7 # bp-rp mag

# Plotting the GCNS magnitudes and the on magnitudes on same plot
fig = plt.figure(figsize = (10,20))

plt.scatter(x_mag_GCNS, gmag_GCNS, c=x_mag_GCNS, cmap=cm.RdBu, s=10) # colourmap using cmap for red-blue
plt.scatter(x_mag_on, gmag_on, color = 'black', marker = '*')
plt.scatter(x_sun, y_sun, color = 'orange', marker = '*')
plt.legend(['All Stars','Stars On Ellipsoid', 'Sun' ])
plt.xlabel('BPmag-RPmag')
plt.ylabel('Gmag')

plt.show()