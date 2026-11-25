# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 15:48:24 2025

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



planet_data = pd.read_csv(r"C:\Users\meg16\Downloads\Data\Data\habitableworldscatalogue.csv")


GB_long = 79.8398
GB_lat = 38.4331
GB_h = 825


spect_filter =  planet_data[planet_data['S_TYPE'].astype(str).str.startswith(('F', 'K', 'G'))]

# Filter for stars visible by the Green Bank Telescope
visible_planets = spect_filter[(spect_filter['P_DISTANCE']*3.26156<=100000) & (spect_filter['S_DEC'] >= -46) & (spect_filter['S_DEC'] <= 90) & (spect_filter['S_RA'] >= 0) & (spect_filter['S_RA'] <=24)]

ra_radians = np.radians(visible_planets['S_RA'] * 15)  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_planets['S_DEC'])


visible_ra= visible_planets['S_RA']

visible_dec = visible_planets['S_DEC']
planet_dist = visible_planets['P_DISTANCE']*3.26156


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


visible_planets['P_DISTANCE'] = visible_planets['P_DISTANCE']*3.26156  # changing dataset so in light years instead of parsecs

target_year = 2001
range_years = 5

# Years since the event
catalog_year = 1987  
years_since = target_year - catalog_year

# Filter the visible_stars DataFrame based on the correct range
visible_planets = visible_planets[(visible_planets['schedule']>= target_year - range_years) & 
                              (visible_planets['schedule'] <= target_year + range_years)]

# converting all visible stars to correct coordinates to then plot
delta_vs_ra = np.radians(visible_planets['S_RA']*15) * years_since  # Convert arcseconds to degrees
delta_vs_dec = np.radians(visible_planets['S_DEC']) * years_since   # Convert arcseconds to degrees


delta_x = visible_planets['P_DISTANCE'] * years_since * np.cos(delta_vs_dec) * np.cos(delta_vs_ra)
delta_y = visible_planets['P_DISTANCE'] * years_since * np.cos(delta_vs_dec) * np.sin(delta_vs_ra)
delta_z = visible_planets['P_DISTANCE'] * years_since * np.sin(delta_vs_dec)


# Recalculate the RA and Dec after applying proper motion
new_vs_ra =visible_planets['S_RA'] + delta_vs_ra
new_vs_dec = visible_planets['S_DEC'] + delta_vs_dec

# Recalculate the Cartesian coordinates
new_vs_ra_radians = np.radians(new_vs_ra * 15)  # Convert from hours to degrees, then to radians
new_vs_dec_radians = np.radians(new_vs_dec)

new_planet_x = visible_planets['P_DISTANCE'] * np.cos(new_vs_dec_radians) * np.cos(new_vs_ra_radians)
new_planet_y = visible_planets['P_DISTANCE'] * np.cos(new_vs_dec_radians) * np.sin(new_vs_ra_radians)
new_planet_z = visible_planets['P_DISTANCE'] * np.sin(new_vs_dec_radians)



target_schedules = target_year
target_schedules_min =  target_year - range_years
target_schedules_max =  target_year + range_years


# Coordinates of the foci
earth_coords = np.array([0, 0, 0])  # Earth's location at origin
event_coords = np.array([x_sn, y_sn, z_sn])

event_dist_2d = np.sqrt(x_sn**2+y_sn**2)


# Assume we have the proper motion in RA (in arcseconds/year) and Dec (in arcseconds/year)
# Also, we have the radial velocity (in km/s)

# For simplicity, let's assume the reference year is the catalog year of the star.
# The star's position is adjusted based on proper motion and radial velocity for the target year.

# Setting up for ellipse
c_sch  = visible_planets.loc[round(visible_planets['schedule']) == target_year]
c_sch1 = visible_planets.loc[round(visible_planets['schedule']) == round(visible_planets['schedule']).max()]
c_sch2 = visible_planets.loc[round(visible_planets['schedule']) == round(visible_planets['schedule']).min()]


# Adjust the RA and Dec for proper motion (assuming pm_ra and pm_dec are in arcseconds/year)
delta_ra2 = np.radians(float(c_sch2['S_RA'].values[0])*15) * years_since  # Convert arcseconds to degrees
delta_dec2 = np.radians(float(c_sch2['S_DEC'].values[0])) * years_since  # Convert arcseconds to degrees

delta_x2 = float(c_sch2['P_DISTANCE'].values[0]) * years_since * np.cos(delta_dec2) * np.cos(delta_ra2)
delta_y2 = float(c_sch2['P_DISTANCE'].values[0]) * years_since * np.cos(delta_dec2) * np.sin(delta_ra2)
delta_z2 = float(c_sch2['P_DISTANCE'].values[0]) * years_since * np.sin(delta_dec2)


# Recalculate the RA and Dec after applying proper motion
new_ra2 = float(c_sch2['S_RA'].values[0]) + delta_ra2
new_dec2 = float(c_sch2['S_DEC'].values[0]) + delta_dec2

# Recalculate the Cartesian coordinates
new_ra_radians2 = np.radians(new_ra2 * 15)  # Convert from hours to degrees, then to radians
new_dec_radians2 = np.radians(new_dec2)

new_planet_x2 = float(c_sch2['P_DISTANCE'].values[0]) * np.cos(new_dec_radians2) * np.cos(new_ra_radians2)
new_planet_y2 = float(c_sch2['P_DISTANCE'].values[0]) * np.cos(new_dec_radians2) * np.sin(new_ra_radians2)
new_planet_z2 = float(c_sch2['P_DISTANCE'].values[0]) * np.sin(new_dec_radians2)


# Calculate distances from the event (e.g., SN1987A)
event_2_starx2 = x_sn - new_planet_x2
event_2_stary2 = y_sn - new_planet_y2
event_2_starz2 = z_sn - new_planet_z2
star_2_earth2 = float(c_sch2['P_DISTANCE'].values[0])

event_2_star2 = np.sqrt((new_planet_x2 - x_sn)**2 + (new_planet_y2 - y_sn)**2)
# Example focii and sum-distance
sn_foci_x,sn_foci_y  = x_sn, y_sn
e_foci_x, e_foci_y = 0, 0

c2 = event_2_star2+star_2_earth2

# Compute ellipse parameters
a2 = c2 / 2                                # Semimajor axis
x0 = (sn_foci_x + e_foci_x) / 2                       # Center x-value
y0 = (sn_foci_y + e_foci_y) / 2                       # Center y-value
f = np.sqrt((sn_foci_x - x0)**2 + (sn_foci_y - y0)**2) # Distance from center to focus
b2 = np.sqrt(a2**2 - f**2)                 # Semiminor axis
phi = np.arctan2((e_foci_y - sn_foci_y), (e_foci_x - sn_foci_x))   # Angle betw major axis and x-axis

# Parametric plot in t
resolution = 10000
t = np.linspace(0, 2*np.pi, resolution)
x2 = x0 + a2 * np.cos(t) * np.cos(phi) - b2 * np.sin(t) * np.sin(phi)
y2 = y0 + a2 * np.cos(t) * np.sin(phi) + b2 * np.sin(t) * np.cos(phi)

# Plot ellipse

plt.plot(x2,y2, color = 'orange')

plt.scatter(new_planet_x, new_planet_y, color = 'purple', marker = '*')
# Show focii
plt.plot(sn_foci_x, sn_foci_y, 'bo')
plt.plot(e_foci_x, e_foci_y, 'go')

plt.axis('equal')
#plt.legend([2025, 'stars', 'Event', 'Earth' ], fontsize = '20')
plt.xlabel('Lightyears', fontsize = '20')
plt.ylabel('Lightyears', fontsize = '20')
plt.xticks(fontsize = '20')
plt.yticks(fontsize = '20')
plt.show()

#%%

planet_dist_new = visible_planets['S_DISTANCE']*3.26156
planet_ra_new = visible_planets['S_RA']*15
planet_dec_new = visible_planets['S_DEC']

dist_max_error = planet_dist_new+(visible_planets['S_DISTANCE_ERROR_MAX']*3.26156)
dist_min_error = planet_dist_new+(visible_planets['S_DISTANCE_ERROR_MIN']*3.26156)


x_sn = dist_sn * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
y_sn = dist_sn  * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
z_sn = dist_sn *np.sin(np.radians(planet_dec_new))

visible_planets_x = planet_dist_new * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
visible_planets_y = planet_dist_new * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
visible_planets_z = planet_dist_new * np.sin(np.radians(planet_dec_new))


max_planets_x = dist_max_error * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
max_planets_y = dist_max_error * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
max_planets_z = dist_max_error * np.sin(np.radians(planet_dec_new))

min_planets_x = dist_min_error * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
min_planets_y = dist_min_error * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
min_planets_z = dist_min_error * np.sin(np.radians(planet_dec_new))

max_sn = dist_sn+(1.2e3*3.26156)
max_x_sn = max_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
max_y_sn = max_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
max_z_sn = max_sn * np.sin(np.radians(dec_sn))

max_d2x = np.sqrt(max_x_sn**2 + max_planets_x**2)
max_d2y = np.sqrt(max_y_sn**2 + max_planets_y**2)
max_d2z = np.sqrt(max_z_sn**2 + max_planets_z**2)

max_d2 = np.sqrt(max_d2x**2+max_d2y**2+max_d2z**2)

min_sn = dist_sn-(1.2e3*3.26156)

min_x_sn = min_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
min_y_sn = min_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
min_z_sn = min_sn * np.sin(np.radians(dec_sn))

min_d2x = np.sqrt(min_x_sn**2 + min_planets_x**2)
min_d2y = np.sqrt(min_y_sn**2 + min_planets_y**2)
min_d2z = np.sqrt(min_z_sn**2 + min_planets_z**2)

min_d2 = np.sqrt(min_d2x**2+min_d2y**2+min_d2z**2)

max_event_dist = np.sqrt(max_x_sn**2+max_y_sn**2+max_z_sn**2)
event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)
min_event_dist = np.sqrt(min_x_sn**2+min_y_sn**2+min_z_sn**2)


d2x = x_sn-visible_planets_x
d2y = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)

earth = 0

max_d3 = np.sqrt(dist_max_error**2+earth**2)
d3 = planet_dist_new-earth
min_d3 = np.sqrt(dist_min_error**2+earth**2)

max_d1 = np.sqrt(max_event_dist**2+earth**2)
d1 = event_dist-earth
min_d1 = np.sqrt(min_event_dist**2+earth**2) 

max_d2_d3 = np.sqrt(max_d2**2+max_d3**2)
min_d2_d3 = np.sqrt(min_d2**2+min_d3**2)


signal_time = (d2+d3)-d1
min_signal_time = np.sqrt(min_d2_d3**2+min_d1**2)#/signal_time
max_signal_time = np.sqrt(max_d2_d3**2+max_d1**2)#/signal_time

min_schedule = 1987+min_signal_time
max_schedule = 1987+max_signal_time
from scipy.interpolate import interp1d

# Create interpolation functions with sorted data
interp_min_func = interp1d(min_signal_time, dist_min_error, kind='linear', fill_value="extrapolate")
interp_max_func = interp1d(max_signal_time, dist_max_error, kind='linear', fill_value="extrapolate")
interp_func = interp1d(signal_time, planet_dist_new, kind='linear', fill_value="extrapolate")

# Generate new interpolated points
interp_time_range = np.linspace(signal_time.min(), signal_time.max(), 100)

interp_dist = interp_func(interp_time_range)

interp_time_range_min = np.linspace(min_signal_time.min(), min_signal_time.max(), 100)
interp_min_dist = interp_min_func(interp_time_range_min)

interp_time_range_max = np.linspace(max_signal_time.min(), max_signal_time.max(), 100)
interp_max_dist = interp_max_func(interp_time_range_max)


# # Create subplots~
fig, ax = plt.subplots(1, 3, figsize=(12, 5))

# Min error plot
ax[0].plot(interp_time_range_min, interp_min_dist,  label='Min Error')
ax[0].set_xlabel('Min Error Time')
ax[0].set_ylabel('Min Error Distance (ly)')
ax[0].legend()
ax[0].grid(True)

# Max error plot
ax[1].plot(interp_time_range_max, interp_max_dist, label='Max Error', color='red')
ax[1].set_xlabel('Max Error Time')
ax[1].set_ylabel('Max Error Distance (ly)')
ax[1].legend()
ax[1].grid(True)

ax[2].plot(interp_time_range, interp_dist, label = 'Original', color='red')
ax[2].set_xlabel('Time')
ax[2].set_ylabel('Distance (ly)')
ax[2].legend()
ax[2].grid(True)
plt.tight_layout()
plt.show()

# plt.plot(interp_time_range_min, interp_min_dist,  label='Min Error')
# plt.plot(interp_time_range_max, interp_max_dist, label='Max Error', color='red')
# plt.plot(interp_time_range, interp_dist, color='red')




schedule = 1987+signal_time

visible_planets['schedule'] = 1987+round(signal_time)