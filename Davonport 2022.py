# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 14:34:56 2025

@author: meg16
"""

from astropy.io.votable import parse
import pandas as pd
#matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import setisynch.gaia as ss
import setisynch.hygdata_v3 as hs
 
# Parse the VOTable
votable = parse(r"C:\Users\meg16\Downloads\Data (1)\Data\GCNS.vot") 
table = votable.get_first_table().to_table()
 
# Convert to Pandas DataFrame
stars_gaia = table.to_pandas()
 
# Extracting the ra, dec, distance and coordinate values as arrays. 
ra_hyg = np.asanyarray(stars_gaia['RA_ICRS']) # in degrees
dec_hyg = np.asanyarray(stars_gaia['DE_ICRS'])
dist_hyg = np.asanyarray(stars_gaia['Dist50']) # in kiloparsecs 
x_GCNS = np.asanyarray(stars_gaia['xcoord50']) # in parsecs 
y_GCNS = np.asanyarray(stars_gaia['ycoord50']) # in parsecs 
z_GCNS = np.asanyarray(stars_gaia['zcoord50']) # in parsecs 
print(f"Number of rows in the CSV: {len(stars_gaia)}")
 
#Checking the first few coloumns
print(stars_gaia.head())

#Creating event for SN1987A

SN1987A = ss.Event(ra=ss.hms2deg(5, 35, 28.03),
                      dec=ss.deg_arcmin_arcsec2deg(-69, 17, 11.79),
                      dist=ss.pc2ly(51.4e3),
                    event_year=ss.date2years('1987-02-24'),
                      event_name='SN1987A')

# Getting SN1987A ra, dec, and distance

ra_sn = SN1987A.ra
dec_sn = SN1987A.dec
dist_sn = SN1987A.dist

# Calculating the coordinates of the event and the distance of the event
x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn *np.sin(np.radians(dec_sn))
event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

# Convert star positions to light-years
x_GCNS_ly = x_GCNS * 3.26  # Convert parsecs to light-years
y_GCNS_ly = y_GCNS * 3.26
z_GCNS_ly = z_GCNS * 3.26

# Calculate the distance between each star and the supernova
d2x = x_sn - x_GCNS_ly
d2y = y_sn - y_GCNS_ly
d2z = z_sn - z_GCNS_ly

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)

# Creating a threshold to only plot the stars 0.1light years from the event
stars_dist_01 = d2 <= 0.1

stars_near_sn = stars_gaia[stars_dist_01]

# Applying the threshold to the stars coordinates
stars_near_sn_x =  np.asanyarray(stars_near_sn['xcoord50']*3.26) # in lightyears

stars_near_sn_y =  np.asanyarray(stars_near_sn['ycoord50']*3.26) # in lightyears

stars_near_sn_z =  np.asanyarray(stars_near_sn['zcoord50']*3.26) # in lightyears

# Plotting the stars x and y coordinates
fig, ax = plt.subplots(figsize=(8, 8))
#plt.scatter(x_sn, y_sn)
plt.scatter(stars_near_sn_x, stars_near_sn_y, color = 'red')
plt.scatter(y_GCNS_ly, z_GCNS_ly, s=1, color='black', label="Stars")  # Small dots for stars
plt.scatter(0, 0, color='blue', marker='x', s=100, label="Earth")  # Earth at (0,0) with blue cross


