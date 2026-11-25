# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 12:03:35 2025

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



#spect_filter =  planet_data[planet_data['S_TYPE'].astype(str).str.startswith(('F', 'K', 'G'))]

# Filter for stars visible by the Green Bank Telescope
visible_planets = gaia[(gaia['Dist50']*3.26156e3<=100000) & (gaia['DE_ICRS'] >= -46) & (gaia['DE_ICRS'] <= 90) & (gaia['RA_ICRS'] >= 0) & (gaia['RA_ICRS'] <=24)]

ra_radians = np.radians(visible_planets['RA_ICRS'] * 15)  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_planets['DE_ICRS'])


planet_dist = visible_planets['Dist50']*3.26156e3


visible_planets_x = planet_dist * np.cos(dec_radians) * np.cos(ra_radians)
visible_planets_y = planet_dist  * np.cos(dec_radians) * np.sin(ra_radians)
visible_planets_z = planet_dist *np.sin(dec_radians)


 
 #Creating event for SN1987A

SN1987A = ss.Event(ra=ss.hms2deg(5, 35, 28.03),
                      dec=ss.deg_arcmin_arcsec2deg(-69, 17, 11.79),
                      dist=ss.pc2ly(51.4e3),
                    event_year=ss.date2years('1987-02-24'),
                      event_name='SN1987A')

# Getting SN1987A ra, dec, and distance using the named_star object:

ra_sn = SN1987A.ra
dec_sn = SN1987A.dec
dist_sn = SN1987A.dist


x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn *np.sin(np.radians(dec_sn))


event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

d2x = x_sn-visible_planets_x
d2y  = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)

earth = 0

d3= planet_dist-earth

d1 = event_dist-earth

signal_time = (d2+d3)-d1

schedule = 1987+signal_time

visible_planets['schedule'] = 1987+round(signal_time)


#%%

planet_dist_new = planet_dist
planet_ra_new = visible_planets['RA_ICRS']*15
planet_dec_new = visible_planets['DE_ICRS']

dist_max_error = (visible_planets['Dist84']-visible_planets['Dist50'])/visible_planets['Dist50']
dist_min_error = (visible_planets['Dist50']-visible_planets['Dist16'])/visible_planets['Dist50']


# Supernova (sn) position (x_sn, y_sn, z_sn)
x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn * np.sin(np.radians(dec_sn))

# Planet position (x, y, z)
visible_planets_x = planet_dist_new * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
visible_planets_y = planet_dist_new * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
visible_planets_z = planet_dist_new * np.sin(np.radians(planet_dec_new))

# Calculate the distance between the supernova and the visible planets (d2)
d2x = x_sn - visible_planets_x
d2y = y_sn - visible_planets_y
d2z = z_sn - visible_planets_z

# Compute the total distance (d2) between the supernova and planets
d2 = np.sqrt(d2x**2 + d2y**2 + d2z**2)


# Maximum and minimum errors for the planet's x, y, z
max_x_sn = (dist_sn + 1.2e3 * 3.26156) * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
max_y_sn = (dist_sn + 1.2e3 * 3.26156) * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
max_z_sn = (dist_sn + 1.2e3 * 3.26156) * np.sin(np.radians(dec_sn))

min_x_sn = (dist_sn - 1.2e3 * 3.26156) * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
min_y_sn = (dist_sn - 1.2e3 * 3.26156) * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
min_z_sn = (dist_sn - 1.2e3 * 3.26156) * np.sin(np.radians(dec_sn))

# Calculate the uncertainties in x, y, z
max_d2x = np.sqrt(max_x_sn**2 + visible_planets_x**2)
max_d2y = np.sqrt(max_y_sn**2 + visible_planets_y**2)
max_d2z = np.sqrt(max_z_sn**2 + visible_planets_z**2)

min_d2x = np.sqrt(min_x_sn**2 + visible_planets_x**2)
min_d2y = np.sqrt(min_y_sn**2 + visible_planets_y**2)
min_d2z = np.sqrt(min_z_sn**2 + visible_planets_z**2)

# Compute the uncertainties in d2 components
max_d2 = np.sqrt(max_d2x**2 + max_d2y**2 + max_d2z**2)
min_d2 = np.sqrt(min_d2x**2 + min_d2y**2 + min_d2z**2)

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

# Print results
print(f"Max Signal Time: {max_signal_time}")
print(f"Min Signal Time: {min_signal_time}")

# Display the uncertainties and the final distances
print(f"Max distance: {max_d2}")
print(f"Min distance: {min_d2}")

min_schedule = 1987+min_signal_time
max_schedule = 1987+max_signal_time

# Create subplots~
fig, ax = plt.subplots(1,2, figsize=(12, 5))

# Min error plot
ax[0].plot(min_signal_time, dist_min_error,  label='Min Error')
ax[0].set_xlabel('Min Error Time')
ax[0].set_ylabel('Min Error Distance')
ax[0].legend()


# Max error plot
ax[1].plot(max_signal_time, dist_max_error ,label='Max Error', color='red')
ax[1].set_xlabel('Max Error Time')
ax[1].set_ylabel('Max Error Distance')
ax[1].legend()

plt.tight_layout()
plt.show()

# plt.plot(interp_time_range_min, interp_min_dist,  label='Min Error', color = 'green')
# plt.plot(interp_time_range_max, interp_max_dist, label='Max Error', color='red')
# plt.plot(interp_time_range, interp_dist, color='blue')




schedule = 1987+signal_time

visible_planets['schedule'] = 1987+round(signal_time)