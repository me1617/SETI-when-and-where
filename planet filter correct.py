# -*- coding: utf-8 -*-
"""
Created on Mon May  5 13:03:59 2025

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
spect_filter = planet_data[planet_data['S_TYPE_TEMP'].isin(['F', 'K', 'G'])]


# Filter for stars visible by the Green Bank Telescope
visible_planets = spect_filter[(spect_filter['S_DISTANCE']*3.26156<=100000) & (spect_filter['S_DEC'] >= -46) & (spect_filter['S_DEC'] <= 90) &
                               (spect_filter['S_RA'] >= -360) & (spect_filter['S_RA'] <=360)]

ra_radians = np.radians(visible_planets['S_RA'])  # Convert degrees to radians
dec_radians = np.radians(visible_planets['S_DEC']) # Convert to radians

# Define ra, dec and distance of the planets
planet_dist = visible_planets['S_DISTANCE']*3.26156

#Calculating the x,y,z position of the stars
visible_planets_x = planet_dist * np.cos(dec_radians) * np.cos(ra_radians)
visible_planets_y = planet_dist  * np.cos(dec_radians) * np.sin(ra_radians)
visible_planets_z = planet_dist *np.sin(dec_radians)
 
 #Creating event for SN1987A using the setisynch gaia code

SN1987A = ss.Event(ra=ss.hms2deg(5, 35, 28.03),
                      dec=ss.deg_arcmin_arcsec2deg(-69, 17, 11.79),
                      dist=ss.pc2ly(51.4e3),
                    event_year=ss.date2years('1987-02-24'),
                      event_name='SN1987A')

# Getting SN1987A ra, dec, and distance

ra_sn = SN1987A.ra
dec_sn = SN1987A.dec
dist_sn = SN1987A.dist

 # Calculating coordinates for the event
x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn * np.sin(np.radians(dec_sn))

# Calculating the event distance in 2D
event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

# Calculating the distance between the event and the stars
d2x = x_sn-visible_planets_x
d2y  = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)

earth = 0

# Calculating the distance between the stars and the earth
d3 = planet_dist-earth

# Calculating the distande between the event and earth
d1 = event_dist

# Calculating the signal time and the schedule and saving it to the dataset
signal_time = (d2+d3)-d1
schedule = 1987+(signal_time) # Putting all the schedule values as a whole number to simplify the process for calculating the ellipse and not choosing a specific date for the calculation
visible_planets['schedule'] = 1987+(signal_time)  # saving the schedule as the exact schedule value so when convert to proper date and time it is accurate

visible_planets['signal_time'] = signal_time 


# Choosing the year we want to look at and allowing a range of 10 years on the plot
target_year = 2025
range_years = 10

# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year

#%%
# Filtering the catalogue bed on the target year
visible_planets = visible_planets[(visible_planets['schedule'] >= target_year - range_years) & 
                              (visible_planets['schedule'] <= target_year + range_years)]
# Setting up for ellipse

# Filtering targets into different year ranges
c_sch1 = visible_planets.loc[(visible_planets['schedule'].astype(int)) == (visible_planets['schedule'].astype(int)).max()]
c_sch  = visible_planets.loc[(visible_planets['schedule'].astype(int)) == target_year]
c_sch2 = visible_planets.loc[(visible_planets['schedule'].astype(int)) == (visible_planets['schedule'].astype(int)).min()]

# Applying the time to the ra and dec
ra_year = (c_sch['S_RA']) * years_since / 3600  # Convert arcseconds to degrees
dec_year = c_sch['S_DEC']* years_since / 3600
new_ra = np.radians(c_sch['S_RA'] + ra_year)  
new_dec = np.radians(c_sch['S_DEC'] + dec_year)

# New position of star for ellipse
new_star_xe = c_sch['S_DISTANCE']*3.26156 * np.cos(new_dec) * np.cos(new_ra)
new_star_ye = c_sch['S_DISTANCE']*3.26156 * np.cos(new_dec) * np.sin(new_ra)
new_star_ze = c_sch['S_DISTANCE']*3.26156 * np.sin(new_dec)

# Calculate distances from event to star and from star to event 
event_2_starx = x_sn - float(new_star_xe.values[0])
event_2_stary = y_sn - float(new_star_ye.values[0])
event_2_starz = z_sn - float(new_star_ze.values[0])
event_2_star = np.sqrt((x_sn- float(new_star_xe.values[0]))**2 + (y_sn- float(new_star_ye.values[0]))**2)
star_2_earth = float(c_sch1['S_DISTANCE'].values[0])*3.26156

# Foci
sn_foci_x,sn_foci_y  = x_sn, y_sn
e_foci_x, e_foci_y = 0, 0
event_dist_2d = np.sqrt(x_sn**2+y_sn**2) # event distance in 2d 

# Ellipse parameters
x0 = (sn_foci_x + e_foci_x) / 2                       # Center x-value
y0 = (sn_foci_y + e_foci_y) / 2                       # Center y-value
f = np.sqrt((sn_foci_x - x0)**2 + (sn_foci_y - y0)**2) # Distance from center to focus

c = event_2_star+star_2_earth

a = c / 2                                # Semimajor axis
b = np.sqrt(a**2 - f**2)                 # Semiminor axis
phi = np.arctan2((e_foci_y - sn_foci_y), (e_foci_x - sn_foci_x))   # Angle between major y-axis and x-axis

resolution = 10000
t = np.linspace(0, 2*np.pi, resolution)

# Coordinates for the ellipse
x = x0 + a * np.cos(t) * np.cos(phi) - b * np.sin(t) * np.sin(phi)
y = y0 + a * np.cos(t) * np.sin(phi) + b * np.sin(t) * np.cos(phi)

fig = plt.figure()
# Plot ellipse
plt.plot(x, y, color = 'red')
plt.scatter(new_star_xe, new_star_ye, color = 'purple')
# Plot foci
plt.plot(sn_foci_x, sn_foci_y, 'bo')
plt.plot(e_foci_x, e_foci_y, 'go')
plt.legend([2025, 'Planets', 'Event', 'Earth' ], fontsize = '20')
plt.xlabel('lightyears')
plt.ylabel('lightyears')
plt.axis('equal')
plt.show()
#%%
#converting the decimal year to specific date and time

from datetime import datetime, timedelta
# Function to convert the decimal schedule to proper date and time
def real_date(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    days_in_year = (datetime(year + 1, 1, 1) - start_of_year).days
    date = start_of_year + timedelta(days=remainder * days_in_year)
    return date

# Applying function to targets in 2025
date = (c_sch['schedule']).apply(real_date)

c_sch['date'] = date # Saving in file

# Saving all data of interest to put in a table
final_data = list(zip(c_sch['P_NAME'], c_sch['S_RA'], c_sch['S_DEC'], 
                      c_sch['S_DISTANCE'], c_sch['schedule'], c_sch['date']))
# Column headers
headers = ["Planet", "ra", "dec", "dist (ly)", "schedule", "date"]

# Create the table
fig, ax = plt.subplots(figsize=(20,25))
ax.set_axis_off()  

table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')
plt.show()

# Adjust font size for all cells
table.auto_set_font_size(False)
table.set_fontsize(15)  # Set your desired font size

# Optionally adjust row/column sizes
table.scale(1, 2)  # (col_width, row_height)
plt.tight_layout()
plt.show()
