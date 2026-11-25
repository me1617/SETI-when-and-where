# -*- coding: utf-8 -*-
"""
Created on Sat Apr 12 11:48:16 2025

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

spect_filter =  star_data[star_data['spect'].astype(str).str.startswith(('F', 'K', 'G'))]


# Filter for stars visible by the Green Bank Telescope
visible_stars = spect_filter[(spect_filter['dist']<=100000) & (spect_filter['dec'] >= -46) & (spect_filter['dec'] <= 90) & (spect_filter['ra'] >= 0) & (spect_filter['ra'] <=24)]


ra_radians = np.radians(visible_stars['ra'] * 15)  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_stars['dec'])


visible_stars['ra']= ra_radians

visible_stars['dec'] = dec_radians
stars_dist = visible_stars['dist']*3.26156


visible_stars_x = stars_dist * np.cos(dec_radians) * np.cos(ra_radians)
visible_stars_y = stars_dist  * np.cos(dec_radians) * np.sin(ra_radians)
visible_stars_z = stars_dist *np.sin(dec_radians)


 
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


event_dist = np.sqrt(x_sn**2+y_sn**2)

d2x = x_sn-visible_stars_x
d2y  = y_sn-visible_stars_y
d2z = z_sn-visible_stars_z

d2 = np.sqrt(d2x**2+d2y**2)

earth = 0

d3= stars_dist-earth

d1 = event_dist-earth

signal_time = (d2+d3)-d1

schedule = 1987+signal_time

visible_stars['schedule'] = 1987+(signal_time)

visible_stars['dist'] = visible_stars['dist']*3.26156  # in light years

target_year = 2025
range_years = 10

# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year

# Filtering the catalogue bed on the target year
visible_stars = visible_stars[(visible_stars['schedule'] >= target_year - range_years) & 
                              (visible_stars['schedule'] <= target_year + range_years)]


# Correct coordinates for targets based on the year of observation

# coverting proper motion to correct coordinate
ra_star_year = visible_stars['ra'] * years_since / 3600  # Convert arcseconds to degrees
dec_star_year = visible_stars['dec'] * years_since / 3600  

# # Calculate ra and dec using proepr motion
# new_ra_star =visible_stars['ra'] + ra_star_year
# new_dec_star = visible_stars['dec'] + dec_star_year

new_ra_star = np.radians(ra_star_year*15)  
new_dec_star = np.radians(dec_star_year)


# new star position
new_star_x = visible_stars['dist'] * np.cos(new_dec_star) * np.cos(new_ra_star)
new_star_y = visible_stars['dist'] * np.cos(new_dec_star) * np.sin(new_ra_star)
new_star_z = visible_stars['dist'] * np.sin(new_dec_star)


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
c_sch1 = visible_stars.loc[round(visible_stars['schedule']) == round(visible_stars['schedule']).max()]
c_sch  = visible_stars.loc[round(visible_stars['schedule']) == target_year]
c_sch2 = visible_stars.loc[round(visible_stars['schedule']) == round(visible_stars['schedule']).min()]

# ra and dec in terms of proper motion
ra_year = float(c_sch['ra'].values[0]) * years_since / 3600  # Convert arcseconds to degrees
dec_year = float(c_sch['dec'].values[0]) * years_since / 3600


new_ra = np.radians(ra_year)  
new_dec = np.radians(dec_year)

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
#%%

# ra and dec in terms of proper motion
ra_year1 = float(c_sch1['ra'].values[0]) * years_since / 3600  # Convert arcseconds to degrees
dec_year1 = float(c_sch1['dec'].values[0]) * years_since / 3600

new_ra1 = ra_year1
new_dec1 = dec_year1

# new position of star for ellipse
new_star_xe1 = float(c_sch1['dist'].values[0]) * np.cos(new_dec1) * np.cos(new_ra1)
new_star_ye1 = float(c_sch1['dist'].values[0]) * np.cos(new_dec1) * np.sin(new_ra1)
new_star_ze1 = float(c_sch1['dist'].values[0]) * np.sin(new_dec1)

# Calculate distances from the event
event_2_starx1 = x_sn - new_star_xe1
event_2_stary1 = y_sn - new_star_ye1
event_2_starz1 = z_sn - new_star_ze1
star_2_earth1 = float(c_sch1['dist'].values[0])

event_2_star1 = np.sqrt((new_star_xe1 - x_sn)**2 + (new_star_ye1 - y_sn)**2)

# foci
sn_foci_x,sn_foci_y  = x_sn, y_sn
e_foci_x, e_foci_y = 0, 0


x0 = (sn_foci_x + e_foci_x) / 2                       # Center x-value
y0 = (sn_foci_y + e_foci_y) / 2                       # Center y-value
f = np.sqrt((sn_foci_x - x0)**2 + (sn_foci_y - y0)**2) # Distance from center to focus

c1 = event_2_star1+star_2_earth1

# Compute ellipse parameters
a1 = c1 / 2                                # Semimajor axis
b1 = np.sqrt(a1**2 - f**2)                 # Semiminor axis
phi = np.arctan2((e_foci_y - sn_foci_y), (e_foci_x - sn_foci_x))   # Angle betw major y-axis and x-axis

# Parametric plot in t
resolution = 10000
t = np.linspace(0, 2*np.pi, resolution)
x1 = x0 + a1 * np.cos(t) * np.cos(phi) - b1 * np.sin(t) * np.sin(phi)
y1 = y0 + a1* np.cos(t) * np.sin(phi) + b1 * np.sin(t) * np.cos(phi)
#%%

# ra and dec in terms of proper motion
ra_year2 = float(c_sch2['ra'].values[0]) * years_since / 3600  # Convert arcseconds to degrees
dec_year2 = float(c_sch2['dec'].values[0]) * years_since / 3600

new_ra2 = float(c_sch2['ra'].values[0]) + ra_year2
new_dec2 = float(c_sch2['dec'].values[0]) + dec_year2

new_ra2 = np.radians(new_ra2 * 15)  
new_dec2 = np.radians(new_dec2)
new_ra2 = ra_year2
new_dec2 = dec_year2
# new position of star for ellipse
new_star_xe2 = float(c_sch2['dist'].values[0]) * np.cos(new_dec2) * np.cos(new_ra2)
new_star_ye2 = float(c_sch2['dist'].values[0]) * np.cos(new_dec2) * np.sin(new_ra2)
new_star_ze2 = float(c_sch2['dist'].values[0]) * np.sin(new_dec2)

# Calculate distances from the event
event_2_starx2= x_sn - new_star_xe2
event_2_stary2 = y_sn - new_star_ye2
event_2_starz2 = z_sn - new_star_ze2
star_2_earth2 = float(c_sch2['dist'].values[0])

event_2_star2 = np.sqrt((new_star_xe2 - x_sn)**2 + (new_star_ye2 - y_sn)**2)

# foci
sn_foci_x,sn_foci_y  = x_sn, y_sn
e_foci_x, e_foci_y = 0, 0


x0 = (sn_foci_x + e_foci_x) / 2                       # Center x-value
y0 = (sn_foci_y + e_foci_y) / 2                       # Center y-value
f = np.sqrt((sn_foci_x - x0)**2 + (sn_foci_y - y0)**2) # Distance from center to focus

c2 = event_2_star2+star_2_earth2

# Compute ellipse parameters
a2 = c2 / 2                                # Semimajor axis
b2 = np.sqrt(a2**2 - f**2)                 # Semiminor axis
phi = np.arctan2((e_foci_y - sn_foci_y), (e_foci_x - sn_foci_x))   # Angle betw major y-axis and x-axis

# Parametric plot in t
resolution = 10000
t = np.linspace(0, 2*np.pi, resolution)
x2 = x0 + a2 * np.cos(t) * np.cos(phi) - b2 * np.sin(t) * np.sin(phi)
y2 = y0 + a2 * np.cos(t) * np.sin(phi) + b2 * np.sin(t) * np.cos(phi)


# Plot ellipse
plt.plot(x, y, color = 'red')
plt.plot(x1, y1, color = 'orange')
plt.plot(x2,y2, color = 'blue')
plt.scatter(new_star_x, new_star_y, color = 'purple')
# Show focii
plt.plot(sn_foci_x, sn_foci_y, 'bo')
plt.plot(e_foci_x, e_foci_y, 'go')
plt.legend([2025, 2035, 2015, 'stars', 'Event', 'Earth' ], fontsize = '20')
plt.xlabel('lightyears')
plt.ylabel('lightyears')
plt.axis('equal')
plt.show()
#%%

import numpy as np
import matplotlib.pyplot as plt

# Define the target year
target_year = 2035

# Extract data for the specific years
c_sch2 = visible_stars.loc[round(visible_stars['schedule']) == round(visible_stars['schedule']).min()]
c_sch  = visible_stars.loc[round(visible_stars['schedule']) == target_year]
c_sch1 = visible_stars.loc[round(visible_stars['schedule']) == round(visible_stars['schedule']).max()]

# Function to convert RA and Dec to Cartesian coordinates
def ra_dec_to_cartesian(ra_deg, dec_deg, distance):

    x = distance * np.cos(dec_deg) * np.cos(ra_deg)
    y = distance * np.cos(dec_deg) * np.sin(ra_deg)
    z = distance * np.sin(dec_deg)
    return x, y, z

# Function to compute ellipse parameters and coordinates
def compute_ellipse(star_data, x_sn, y_sn):
    ra = float(star_data['ra'].values[0])
    dec = float(star_data['dec'].values[0])
    dist = float(star_data['dist'].values[0])
    
    # Convert to Cartesian coordinates
    x_star, y_star, z_star = ra_dec_to_cartesian(ra, dec, dist)
    
    # Calculate distances
    event_to_star = np.sqrt((x_star - x_sn)**2 + (y_star - y_sn)**2)
    star_to_earth = dist
    c = event_to_star + star_to_earth  # Total distance (2a)
    
    # Foci positions
    sn_foci_x, sn_foci_y = x_sn, y_sn
    e_foci_x, e_foci_y = 0, 0
    
    # Center of the ellipse
    x0 = (sn_foci_x + e_foci_x) / 2
    y0 = (sn_foci_y + e_foci_y) / 2
    
    # Distance from center to focus
    f = np.sqrt((sn_foci_x - x0)**2 + (sn_foci_y - y0)**2)
    
    # Semi-major and semi-minor axes
    a = c / 2
    b = np.sqrt(a**2 - f**2)
    
    # Rotation angle
    phi = np.arctan2((e_foci_y - sn_foci_y), (e_foci_x - sn_foci_x))
    
    # Parametric equations for the ellipse
    resolution = 10000
    t = np.linspace(0, 2*np.pi, resolution)
    x = x0 + a * np.cos(t) * np.cos(phi) - b * np.sin(t) * np.sin(phi)
    y = y0 + a * np.cos(t) * np.sin(phi) + b * np.sin(t) * np.cos(phi)
    
    return x, y, x_star, y_star


x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn *np.sin(np.radians(dec_sn))


# Compute ellipses for each observation year
x2, y2, x_star2, y_star2 = compute_ellipse(c_sch2, x_sn, y_sn)
x, y, x_star, y_star = compute_ellipse(c_sch, x_sn, y_sn)
x1, y1, x_star1, y_star1 = compute_ellipse(c_sch1, x_sn, y_sn)

# Plot the ellipses and key points
plt.figure(figsize=(10, 8))
plt.plot(x, y, color='red', label=f'{target_year}')
plt.plot(x1, y1, color='orange', label='Latest Year')
plt.plot(x2, y2, color='orange', label='Earliest Year')
plt.scatter(x_star, y_star, color='purple', label='Star Position')
plt.plot(x_sn, y_sn, 'bo', label='Supernova Event')
plt.plot(0, 0, 'go', label='Earth')
plt.xlabel('Lightyears')
plt.ylabel('Lightyears')
#plt.axis('equal')
plt.legend(fontsize=12)
plt.title('Ellipses Representing Light Travel Paths')
plt.grid(True)
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