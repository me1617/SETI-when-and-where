# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 17:16:11 2025

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


# Load the planet catalog
planet_data = pd.read_csv(r"C:\Users\meg16\Downloads\Data\Data\habitableworldscatalogue.csv")


#Filter the planet data based on only on the habitable planets
spect_filter =  planet_data[planet_data['S_TYPE'].astype(str).str.startswith(('F', 'K', 'G'))]

# Filter for stars visible by the Green Bank Telescope
visible_planets = spect_filter[(spect_filter['S_DISTANCE']*3.26156<=100000) & (spect_filter['S_DEC'] >= -46) & (spect_filter['S_DEC'] <= 90) & (spect_filter['S_RA'] >= 0) & (spect_filter['S_RA'] <=24)]

ra_radians = np.radians(visible_planets['S_RA'] * 15)  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_planets['S_DEC']) # Convert to radians

# Define ra, dec and distance of the planets
visible_ra= visible_planets['S_RA']
visible_dec = visible_planets['S_DEC']
planet_dist = visible_planets['S_DISTANCE']*3.26156

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

visible_planets['schedule'] = 1987+(signal_time)


visible_planets['S_DISTANCE'] = visible_planets['S_DISTANCE']*3.26156  # changing dataset so in light years instead of parsecs

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


delta_x = visible_planets['S_DISTANCE'] * years_since * np.cos(delta_vs_dec) * np.cos(delta_vs_ra)
delta_y = visible_planets['S_DISTANCE'] * years_since * np.cos(delta_vs_dec) * np.sin(delta_vs_ra)
delta_z = visible_planets['S_DISTANCE'] * years_since * np.sin(delta_vs_dec)


# Recalculate the RA and Dec after applying proper motion
new_vs_ra =visible_planets['S_RA'] + delta_vs_ra
new_vs_dec = visible_planets['S_DEC'] + delta_vs_dec

# Recalculate the Cartesian coordinates
new_vs_ra_radians = np.radians(new_vs_ra * 15)  # Convert from hours to degrees, then to radians
new_vs_dec_radians = np.radians(new_vs_dec)

new_planet_x = visible_planets['S_DISTANCE'] * np.cos(new_vs_dec_radians) * np.cos(new_vs_ra_radians)
new_planet_y = visible_planets['S_DISTANCE'] * np.cos(new_vs_dec_radians) * np.sin(new_vs_ra_radians)
new_planet_z = visible_planets['S_DISTANCE'] * np.sin(new_vs_dec_radians)



target_schedules = target_year
target_schedules_min =  target_year - range_years
target_schedules_max =  target_year + range_years


# Coordinates of the foci
earth_coords = np.array([0, 0, 0])  # Earth's location at origin
event_coords = np.array([x_sn, y_sn, z_sn])

event_dist_2d = np.sqrt(x_sn**2+y_sn**2)


# Setting up for ellipse
c_sch  = visible_planets.loc[round(visible_planets['schedule']) == target_year]
c_sch1 = visible_planets.loc[round(visible_planets['schedule']) == round(visible_planets['schedule']).max()]
c_sch2 = visible_planets.loc[round(visible_planets['schedule']) == round(visible_planets['schedule']).min()]


# Adjust the RA and Dec for proper motion (assuming pm_ra and pm_dec are in arcseconds/year)
delta_ra2 = np.radians(float(c_sch1['S_RA'].values[0])*15) * years_since  # Convert arcseconds to degrees
delta_dec2 = np.radians(float(c_sch1['S_DEC'].values[0])) * years_since  # Convert arcseconds to degrees

delta_x2 = float(c_sch1['S_DISTANCE'].values[0]) * years_since * np.cos(delta_dec2) * np.cos(delta_ra2)
delta_y2 = float(c_sch1['S_DISTANCE'].values[0]) * years_since * np.cos(delta_dec2) * np.sin(delta_ra2)
delta_z2 = float(c_sch1['S_DISTANCE'].values[0]) * years_since * np.sin(delta_dec2)


# Recalculate the RA and Dec after applying proper motion
new_ra2 = float(c_sch1['S_RA'].values[0]) + delta_ra2
new_dec2 = float(c_sch1['S_DEC'].values[0]) + delta_dec2

# Recalculate the Cartesian coordinates
new_ra_radians2 = np.radians(new_ra2 * 15)  # Convert from hours to degrees, then to radians
new_dec_radians2 = np.radians(new_dec2)

new_planet_x2 = float(c_sch1['S_DISTANCE'].values[0]) * np.cos(new_dec_radians2) * np.cos(new_ra_radians2)
new_planet_y2 = float(c_sch1['S_DISTANCE'].values[0]) * np.cos(new_dec_radians2) * np.sin(new_ra_radians2)
new_planet_z2 = float(c_sch1['S_DISTANCE'].values[0]) * np.sin(new_dec_radians2)


# Calculate distances from the event (e.g., SN1987A)
event_2_starx2 = x_sn - new_planet_x2
event_2_stary2 = y_sn - new_planet_y2
event_2_starz2 = z_sn - new_planet_z2
star_2_earth2 = float(c_sch1['S_DISTANCE'].values[0])

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

plt.scatter(new_planet_x2, new_planet_y2, color = 'purple', marker = '*')
# Show focii
plt.plot(sn_foci_x, sn_foci_y, 'bo')
plt.plot(e_foci_x, e_foci_y, 'go')

plt.axis('equal')
plt.legend([2030, 'stars', 'Event', 'Earth' ], fontsize = '20')
plt.xlabel('Lightyears', fontsize = '20')
plt.ylabel('Lightyears', fontsize = '20')
plt.xticks(fontsize = '20')
plt.yticks(fontsize = '20')
plt.show()


#converting the decimal year to specific date

from datetime import datetime, timedelta

def decimal_year_to_date(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    days_in_year = (datetime(year + 1, 1, 1) - start_of_year).days
    date = start_of_year + timedelta(days=remainder * days_in_year)
    return date
#%%

date = (c_sch1['schedule']).apply(decimal_year_to_date)

c_sch1['date'] = date

final_data = list(zip(c_sch1['S_TYPE'], c_sch1['S_RA'], c_sch1['S_DEC'], 
                      c_sch1['S_DISTANCE'], c_sch1['schedule'], c_sch1['date']))
# Column headers
headers = ["S_TYPE", "S_RA", "S_DEC", "S_DISTANCE (ly)", "schedule", "date"]

# Creating the table
fig, ax = plt.subplots(figsize=(5, 3))
ax.set_axis_off() 

table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')

# Set font size
table.auto_set_font_size(False)  
table.set_fontsize(10)  # Change the size as needed



# Show the figure
plt.show()




