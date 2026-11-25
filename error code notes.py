# -*- coding: utf-8 -*-
"""
Created on Thu May  1 09:42:59 2025

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


print(os.getcwd())
sys.path.append("/Users/meg16/Downloads/Data (2)/Data/Ellipsoid_Functions.py")

# Parse the VOTable
votable = parse(r"C:\Users\meg16\Downloads\Data (2)\Data\GCNS_morecol.vot")  # Update with your actual file path
table = votable.get_first_table().to_table()
# Convert to Pandas DataFrame
gaia = table.to_pandas()
print(gaia.columns)


# Filter for stars visible by the Green Bank Telescope
visible_planets = gaia[(gaia['Dist50']*3.26156e3<=100000) & (gaia['DE_ICRS'] >= -46) & (gaia['DE_ICRS'] <= 90) & (gaia['RA_ICRS'] >= 0) & (gaia['RA_ICRS'] <=24)]

# Load the RA, dec and distance values
ra_radians = np.radians(visible_planets['RA_ICRS'] * 15)  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_planets['DE_ICRS'])

planet_dist = visible_planets['Dist50']*3.26156

# Calculate the coordinates of the planets
visible_planets_x = planet_dist * np.cos(dec_radians) * np.cos(ra_radians)
visible_planets_y = planet_dist  * np.cos(dec_radians) * np.sin(ra_radians)
visible_planets_z = planet_dist *np.sin(dec_radians)
 
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

# Calculate the coordinates of the event and the distance
x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn *np.sin(np.radians(dec_sn))
event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

# Calculate the coordinates and distance between the event and planets
d2x = x_sn-visible_planets_x
d2y  = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z
d2 = np.sqrt(d2x**2+d2y**2)

earth = 0

#Calculate the distance between the planets and earth
d3= planet_dist-earth

# Calculate the distance between the event and earth
d1 = event_dist-earth

#Calculate the difference in time that the planets saw the event comapred to us and calculating the schedule 
signal_time = (d2+d3)-d1

schedule = 1987+signal_time

visible_planets['schedule'] = 1987+round(signal_time) # Save schedule to the file


#%%

#Extracting the planet distances, ra and dec and calculate the corresponding errors
planet_dist_new = planet_dist
planet_ra_new = visible_planets['RA_ICRS']*15
planet_dec_new = visible_planets['DE_ICRS']


dist_max_error = (visible_planets['Dist84']-visible_planets['Dist50'])/visible_planets['Dist50']
dist_min_error = (visible_planets['Dist50']-visible_planets['Dist16'])/visible_planets['Dist50']

# Calculating the coordinates for the event 
x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn *np.sin(np.radians(dec_sn))

# Calculate the new coordinates of the planets and their corresponding maximum and minimum errors
visible_planets_x = planet_dist_new * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
visible_planets_y = planet_dist_new * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
visible_planets_z = planet_dist_new * np.sin(np.radians(planet_dec_new))


max_planets_x = dist_max_error * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
max_planets_y = dist_max_error * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
max_planets_z = dist_max_error * np.sin(np.radians(planet_dec_new))

min_planets_x = dist_min_error * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
min_planets_y = dist_min_error * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
min_planets_z = dist_min_error * np.sin(np.radians(planet_dec_new))

# Calculate the maximum and minimum errors on the event, and distance error between the event and planets
max_sn = dist_sn+(1.2e3*3.26156)
max_x_sn = max_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
max_y_sn = max_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
max_z_sn = max_sn * np.sin(np.radians(dec_sn))

max_d2x = np.sqrt(max_x_sn**2 + max_planets_x**2) # Using propagating error equation
max_d2y = np.sqrt(max_y_sn**2 + max_planets_y**2)
max_d2z = np.sqrt(max_z_sn**2 + max_planets_z**2)

valid_indices1 = ~np.isnan(max_d2x) & ~np.isnan(max_d2y) & ~np.isnan(max_d2z) # Find non-NaN values ensuring no errors
max_d2x = max_d2x[valid_indices1]
max_d2y = max_d2y[valid_indices1]
max_d2z = max_d2z[valid_indices1]

max_d2 = np.sqrt(max_d2x**2+max_d2y**2+max_d2z**2)

min_sn = dist_sn-(1.2e3*3.26156)

min_x_sn = min_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
min_y_sn = min_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
min_z_sn = min_sn * np.sin(np.radians(dec_sn))

min_d2x = np.sqrt(min_x_sn**2 + min_planets_x**2)
min_d2y = np.sqrt(min_y_sn**2 + min_planets_y**2)
min_d2z = np.sqrt(min_z_sn**2 + min_planets_z**2)


valid_indices2 = ~np.isnan(min_d2x) & ~np.isnan(min_d2y) & ~np.isnan(min_d2z) # Find non-NaN values
min_d2x = min_d2x[valid_indices2]
min_d2y = min_d2y[valid_indices2]
min_d2z = min_d2z[valid_indices2]

min_d2 = np.sqrt(min_d2x**2+min_d2y**2+min_d2z**2)

d2x = x_sn - visible_planets_x
d2y = y_sn - visible_planets_y
d2z = z_sn - visible_planets_z

d2 = np.sqrt(d2x**2 + d2y**2 + d2z**2)


# Calculating the maximum and minimum distance for the event
max_event_dist = np.sqrt(max_x_sn**2+max_y_sn**2+max_z_sn**2)
event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)
min_event_dist = np.sqrt(min_x_sn**2+min_y_sn**2+min_z_sn**2)

earth = 0

# Calculating the distance error between the planets and the earth
d3 = planet_dist_new-earth
max_d3 = np.sqrt(dist_max_error**2)
min_d3 = np.sqrt(dist_min_error**2)

# Calculating the distance error between the event and the earth
d1 = event_dist-earth
max_d1 = np.sqrt(max_event_dist**2)
min_d1  = np.sqrt(min_event_dist**2)

# Calculating the signal time error
max_d2_d3 = np.sqrt(max_d2**2+max_d3**2)
min_d2_d3 = np.sqrt(min_d2**2+min_d3**2)
signal_time = (d2+d3)-d1
min_signal_time = ((min_d2 + min_d3) - min_d1) 
max_signal_time = ((max_d2 + max_d3) - max_d1) 
valid_indices5 = ~np.isnan(min_signal_time) & ~np.isnan(max_signal_time) # Find non-NaN values
min_signal_time = min_signal_time[valid_indices5]
max_signal_time = max_signal_time[valid_indices5]


# Calculating the schedule error
min_schedule = 1987+min_signal_time
max_schedule = 1987+max_signal_time

# Using the interpolate library to obtain more points to get the plot smoother and applying it to the signal time and distance errors
from scipy.interpolate import interp1d

# Interpolating
interp_min_func = interp1d(min_signal_time, dist_min_error, kind='linear', fill_value="extrapolate")
interp_max_func = interp1d(max_signal_time, dist_max_error, kind='linear', fill_value="extrapolate")
interp_func = interp1d(signal_time, planet_dist_new, kind='linear', fill_value="extrapolate")

# Generate new interpolated points
interp_time_range = np.linspace(signal_time.min(), signal_time.max(),1000)

interp_dist = interp_func(interp_time_range)

interp_time_range_min = np.linspace(min_signal_time.min(), min_signal_time.max(), 1000)
interp_min_dist = interp_min_func(interp_time_range_min)

interp_time_range_max = np.linspace(max_signal_time.min(), max_signal_time.max(), 1000)
interp_max_dist = interp_max_func(interp_time_range_max)


# Plotting the signal time errors against the distance errors as a subplot
fig, ax = plt.subplots(1, 3, figsize=(12, 5))

# Min error plot
ax[0].plot(interp_time_range_min, interp_min_dist,  label='Min Error')
ax[0].set_xlabel('Min Error Time')
ax[0].set_ylabel('Min Error Distance')
ax[0].legend()
ax[0].grid(True)

# Max error plot
ax[1].plot(interp_time_range_max, interp_max_dist ,label='Max Error', color='red')
ax[1].set_xlabel('Max Error Time')
ax[1].set_ylabel('Max Error Distance')
ax[1].legend()
ax[1].grid(True)

ax[2].plot(interp_time_range, interp_dist ,color='red')
ax[2].set_xlabel('Time')
ax[2].set_ylabel('Distance')
ax[2].legend()
ax[2].grid(True)
plt.tight_layout()
plt.show()

# plt.plot(interp_time_range_min, interp_min_dist,  label='Min Error', color = 'green')
# plt.plot(interp_time_range_max, interp_max_dist, label='Max Error', color='red')
# plt.plot(interp_time_range, interp_dist, color='blue')




schedule = 1987+signal_time

visible_planets['schedule'] = 1987+round(signal_time)