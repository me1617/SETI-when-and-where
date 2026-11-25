# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 15:25:44 2025

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




planet_data = pd.read_csv(r"C:\Users\meg16\Downloads\Data (2)\Data\seen_SN.csv")


# Filter for stars visible by the Green Bank Telescope
visible_planets = planet_data[(planet_data['Dist50']*3.26156e3<=100000) & (planet_data['DE_ICRS'] >= -46) & (planet_data['DE_ICRS'] <= 90) & (planet_data['RA_ICRS'] >= 0) & (planet_data['RA_ICRS'] <=24)]

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
y_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn * np.sin(np.radians(dec_sn))


event_dist = np.sqrt(x_sn**2+y_sn**2)

d2x = x_sn-visible_planets_x
d2y  = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z

d2 = np.sqrt(d2x**2+d2y**2)

earth = 0

d3 = planet_dist-earth

d1 = event_dist-earth

signal_time = (d2+d3)-d1
schedule = 1987+signal_time

visible_planets['schedule'] = 1987+(signal_time)

visible_planets['signal_time'] = signal_time
visible_planets['Dist50'] = visible_planets['Dist50']*3.26156e3 # changing dataset so in light years instead of parsecs
visible_planets['Dist84'] = visible_planets['Dist84']*3.26156e3 
visible_planets['Dist16'] = visible_planets['Dist16']*3.26156e3 

#%%
target_year = np.random.choice(visible_planets['schedule'])
range_years = 10

# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year

# Filtering the catalogue bed on the target year
visible_planets = visible_planets[(visible_planets['schedule'] >= target_year - range_years) & 
                              (visible_planets['schedule'] <= target_year + range_years)]


# Correct coordinates for targets based on the year of observation

# coverting proper motion to correct coordinate
ra_star_year = visible_planets['RA_ICRS'] * years_since / 3600  # Convert arcseconds to degrees
dec_star_year = visible_planets['DE_ICRS'] * years_since / 3600  

# Calculate ra and dec using proepr motion
new_ra_star =visible_planets['RA_ICRS'] + ra_star_year
new_dec_star = visible_planets['DE_ICRS'] + dec_star_year

new_ra_star = np.radians(new_ra_star * 15)  
new_dec_star = np.radians(new_dec_star)

# new star position
new_star_x = visible_planets['Dist50'] * np.cos(new_dec_star) * np.cos(new_ra_star)
new_star_y = visible_planets['Dist50'] * np.cos(new_dec_star) * np.sin(new_ra_star)
new_star_z = visible_planets['Dist50'] * np.sin(new_dec_star)


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
c_sch1 = visible_planets.loc[visible_planets['schedule'] == visible_planets['schedule'].max()]
c_sch  = visible_planets.loc[visible_planets['schedule'] == target_year]
c_sch2 = visible_planets.loc[visible_planets['schedule'] == visible_planets['schedule'].min()]

# ra and dec in terms of proper motion
ra_year = float(c_sch['RA_ICRS'].values[0]) * years_since / 3600  # Convert arcseconds to degrees
dec_year = float(c_sch['DE_ICRS'].values[0]) * years_since / 3600

new_ra = float(c_sch['RA_ICRS'].values[0]) + ra_year
new_dec = float(c_sch['DE_ICRS'].values[0]) + dec_year

new_ra = np.radians(new_ra * 15)  
new_dec = np.radians(new_dec)

# new position of star for ellipse
new_star_xe = float(c_sch['Dist50'].values[0]) * np.cos(new_dec) * np.cos(new_ra)
new_star_ye = float(c_sch['Dist50'].values[0]) * np.cos(new_dec) * np.sin(new_ra)
new_star_ze = float(c_sch['Dist50'].values[0]) * np.sin(new_dec)

# Calculate distances from the event
event_2_starx = x_sn - new_star_xe
event_2_stary = y_sn - new_star_ye
event_2_starz = z_sn - new_star_ze
star_2_earth = float(c_sch['Dist50'].values[0])

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
#converting the decimal year to specific date

from datetime import datetime, timedelta

def decimal_year_to_date(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    days_in_year = (datetime(year + 1, 1, 1) - start_of_year).days
    date = start_of_year + timedelta(days=remainder * days_in_year)
    return date


date = (visible_planets['schedule']).apply(decimal_year_to_date)

c_sch['date'] = date

visible_planets['date'] = date

final_data = list(zip(visible_planets['RA_ICRS'], visible_planets['DE_ICRS'], 
                      visible_planets['Dist50'],  visible_planets['schedule'], visible_planets['date']))
# Column headers
headers = [ "ra", "dec", "dist (ly)","schedule", "date"]
# Create a figure and axis
fig, ax = plt.subplots(figsize=(20,25))
ax.set_axis_off()  


# Create a table
table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')
plt.show()

visible_planets.to_csv('visible_planets.csv', index=False)

#%%
magnitude = visible_planets['Gmag']
signal_time = visible_planets['signal_time']
plt.plot(signal_time, magnitude)

#%%

SN1987A = ss.Event(ra=ss.hms2deg(5, 35, 28.03),
                      dec=ss.deg_arcmin_arcsec2deg(-69, 17, 11.79),
                      dist=ss.pc2ly(51.4e3),
                    event_year=ss.date2years('1987-02-24'),
                      event_name='SN1987A')

df, alerts_scanned = SN1987A.search_alerts(ra_range =  np.max(visible_planets['RA_ICRS']),
                                         dec_range= np.max(visible_planets['DE_ICRS']),
                                         try_synch=False,
                                         scan_type='random',
                                         sample_size=100,
                                         save_df=True,
                                         filename='visible_planets')

#%%

results = []  # Store results

for _, row in visible_planets.iterrows():  # Loop through each row in the DataFrame
    event = ss.Event(
        ra=row['RA_ICRS'],
        dec=row['DE_ICRS'],
        dist=row['Dist50'],
        event_year=row['date'],
        event_name=f"visible_planet_{row.name}"  # Unique event name
    )

    df, alerts_scanned = SN1987A.search_alerts(
        ra_range=row['RA_ICRS'],
        dec_range=row['DE_ICRS'],
        try_synch=False,
        scan_type='random',
        sample_size=100,
        save_df=True,
        filename=f"visible_planet_{row.name}"  # Save each file with unique name
    )

    results.append((df, alerts_scanned))  # Store results for later use

# Now, `results` contains a list of (df, alerts_scanned) for each target




