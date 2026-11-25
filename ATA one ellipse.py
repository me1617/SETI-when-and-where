# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 21:35:11 2025

@author: meg16
"""


import numpy as np
import pandas as pd
import setisynch.gaia as ss
import setisynch.hygdata_v3 as hs
import pyvo as vo
import matplotlib.pyplot as plt
import time
import os
import random



star_data = pd.read_csv(r"C:\Users\meg16\Downloads\Data\Data\athyg_32_hyg_ids.csv")

spect_filter =  star_data[star_data['spect'].astype(str).str.startswith(('G', 'K', 'F'))]

ATA_long = 121.469 # degrees, west
ATA_lat = 40.817  # degrees, north
ATA_h = 1000 # in meters

a = 6378137 # Earth's equitorial plane

b = 6356752.3142 # Earth's polar radius

e_sq = 1-(b**2/a**2) # Earth's eccentricity squared

ATA_N = a/ np.sqrt(1-e_sq*np.sin(ATA_lat))
ATA_x = (ATA_N+ATA_h)*np.cos(ATA_lat)*np.cos(ATA_long)
ATA_y = (ATA_N+ATA_h)*np.cos(ATA_lat)*np.sin(ATA_long)
ATA_z = ((b**2/a**2)*ATA_N+ATA_h)*np.sin(ATA_lat)

# Remove missing data
spect_filter = spect_filter.dropna(subset=['ra', 'dec', 'dist'])

spect_filter_x  = spect_filter['x0']
spect_filter_y  = spect_filter['y0']
spect_filter_z  = spect_filter['z0']

# Filter for stars visible by the Green Bank Telescope


visible_stars_ATA = spect_filter[(spect_filter['dec']*3.26156 >= np.radians(-49.2)) & (spect_filter['dec'] <=np.radians(90)) & (spect_filter['ra'] >= np.radians(0)) & (spect_filter['ra'] <=np.radians(24))]

ra_radians = np.radians(visible_stars_ATA['ra'] * 15)  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_stars_ATA['dec'])
 

visible_stars_x = visible_stars_ATA['x0']
visible_stars_y  = visible_stars_ATA['y0']
visible_stars_z  = visible_stars_ATA['z0']

stars_dist = visible_stars_ATA['dist']*3.26156

#%%
 #Creating event for ETA

ETA = ss.Event(ra=ss.hms2deg(10, 45, 3.591),
                      dec=ss.deg_arcmin_arcsec2deg(-59, 41, 4.26),
                      dist=ss.pc2ly(2300),
                      event_name='ETA')

# Getting ETA ra, dec, and distance using the named_star object:

ra_sn = ETA.ra
dec_sn = ETA.dec
dist_sn = ETA.dist


x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn *np.sin(np.radians(dec_sn))


event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

d2x = x_sn-visible_stars_x
d2y  = y_sn-visible_stars_y
d2z = z_sn-visible_stars_z

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)

earth = 0

d3= stars_dist-earth

d1 = event_dist-earth

signal_time = (d2+d3)-d1

schedule = 1987+signal_time

visible_stars_ATA['schedule'] = 1987+(signal_time)

visible_stars_ATA['dist'] = visible_stars_ATA['dist']*3.26156  # in light years

target_year = 2191
range_years = 10

# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year

# Filtering the catalogue bed on the target year
visible_stars_ATA = visible_stars_ATA[(visible_stars_ATA['schedule'] >= target_year - range_years) & 
                              (visible_stars_ATA['schedule'] <= target_year + range_years)]


# Correct coordinates for targets based on the year of observation

# coverting proper motion to correct coordinate
ra_star_year = visible_stars_ATA['ra'] * years_since / 3600  # Convert arcseconds to degrees
dec_star_year = visible_stars_ATA['dec'] * years_since / 3600  

# Calculate ra and dec using proepr motion
new_ra_star =visible_stars_ATA['ra'] + ra_star_year
new_dec_star = visible_stars_ATA['dec'] + dec_star_year

new_ra_star = np.radians(new_ra_star * 15)  
new_dec_star = np.radians(new_dec_star)

# new star position
new_star_x = visible_stars_ATA['dist'] * np.cos(new_dec_star) * np.cos(new_ra_star)
new_star_y = visible_stars_ATA['dist'] * np.cos(new_dec_star) * np.sin(new_ra_star)
new_star_z = visible_stars_ATA['dist'] * np.sin(new_dec_star)


# starting to create the ellipse
target_schedules = target_year
target_schedules_min =  target_year - range_years
target_schedules_max =  target_year + range_years


# Coordinates of the foci
earth_coords = np.array([0, 0, 0])  # Earth's location at origin
event_coords = np.array([x_sn, y_sn, z_sn])

event_dist_2d = np.sqrt(x_sn**2+y_sn**2) # event distance in 2d 

# Setting up for ellipse

# Correct coordinates for targets based on the year of observation for plotting the ellipse
c_sch = visible_stars_ATA.loc[round(visible_stars_ATA['schedule']) == round(visible_stars_ATA['schedule']).max()]
c_sch1  = visible_stars_ATA.loc[round(visible_stars_ATA['schedule']) == target_year]
c_sch2 = visible_stars_ATA.loc[round(visible_stars_ATA['schedule']) == round(visible_stars_ATA['schedule']).min()]

# ra and dec in terms of proper motion
ra_year = float(c_sch['ra'].values[0]) * years_since / 3600  # Convert arcseconds to degrees
dec_year = float(c_sch['dec'].values[0]) * years_since / 3600

new_ra = float(c_sch['ra'].values[0]) + ra_year
new_dec = float(c_sch['dec'].values[0]) + dec_year

new_ra = np.radians(new_ra * 15)  
new_dec = np.radians(new_dec)

# new position of star for ellipse
new_star_xe = float(c_sch['dist'].values[0]) * np.cos(new_dec) * np.cos(new_ra)
new_star_ye = float(c_sch['dist'].values[0]) * np.cos(new_dec) * np.sin(new_ra)
new_star_ze = float(c_sch['dist'].values[0]) * np.sin(new_dec)

# Calculate distances from the event
event_2_starx = x_sn - new_star_xe
event_2_stary = y_sn - new_star_ye
event_2_starz = z_sn - new_star_ze
star_2_earth = float(c_sch['dist'].values[0])

event_2_star = np.sqrt((new_star_xe - x_sn)**2 + (new_star_ye - y_sn)**2)

# foci
sn_foci_x,sn_foci_y  = x_sn, y_sn
e_foci_x, e_foci_y = 0, 0


x0 = (sn_foci_x + e_foci_x) / 2                       # Center x-value
y0 = (sn_foci_y + e_foci_y) / 2                       # Center y-value
f = np.sqrt((sn_foci_x - x0)**2 + (sn_foci_y - y0)**2) # Distance from center to focus

c = event_2_star+star_2_earth

# Compute ellipse parameters
a = c / 2                                # Semimajor axis
b = np.sqrt(a**2 - f**2)                 # Semiminor axis
phi = np.arctan2((e_foci_y - sn_foci_y), (e_foci_x - sn_foci_x))   # Angle betw major y-axis and x-axis

# Parametric plot in t
resolution = 10000
t = np.linspace(0, 2*np.pi, resolution)
x = x0 + a * np.cos(t) * np.cos(phi) - b * np.sin(t) * np.sin(phi)
y = y0 + a * np.cos(t) * np.sin(phi) + b * np.sin(t) * np.cos(phi)

# Plot ellipse
plt.plot(x, y, color = 'red')

plt.scatter(new_star_x, new_star_y, color = 'purple')
# Show focii
plt.plot(sn_foci_x, sn_foci_y, 'bo')
plt.plot(e_foci_x, e_foci_y, 'go')
plt.legend([target_year, 'stars', 'Event', 'Earth' ], fontsize = '20')
plt.xlabel('lightyears')
plt.ylabel('lightyears')
plt.axis('equal')
plt.show()
#%%
from datetime import datetime, timedelta

def decimal_year_to_date(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    days_in_year = (datetime(year + 1, 1, 1) - start_of_year).days
    date = start_of_year + timedelta(days=remainder * days_in_year)
    return date


date = (c_sch['schedule']).apply(decimal_year_to_date)

c_sch['date'] = date

final_data = list(zip(c_sch['id'], c_sch['ra'], c_sch['dec'], 
                      c_sch['dist'], c_sch['schedule'], c_sch['date']))
# Column headers
headers = ["star_id", "ra", "dec", "dist (ly)", "schedule", "date"]


# Create a figure and axis
fig, ax = plt.subplots(figsize=(20,25))
ax.set_axis_off()  


# Create a table
table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')
plt.show()


