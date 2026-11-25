# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 17:16:11 2025 not

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



visible_planets = pd.read_csv(r"C:\Users\meg16\Downloads\Data\Data\habitableworldscatalogue.csv")


GB_long = 79.8398
GB_lat = 38.4331
GB_h = 825


#spect_filter =  planet_data[planet_data['S_TYPE'].astype(str).str.startswith(('F', 'K', 'G'))]

# Filter for stars visible by the Green Bank Telescope
#visible_planets = planet_data[(planet_data['S_DISTANCE']*3.26156<=100000) & (planet_data['S_DEC'] >= -46) & (planet_data['S_DEC'] <= 90) & (planet_data['S_RA'] >= 0) & (planet_data['S_RA'] <=24)]

ra_radians = np.radians(visible_planets['S_RA'] * 15)  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_planets['S_DEC'])

visible_dec = visible_planets['S_DEC']
planet_dist = visible_planets['S_DISTANCE']*3.26156


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

# Calculating the coordinates for the event
x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn *np.sin(np.radians(dec_sn))


event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

#Calculating the distance between the event and the planets
d2x = x_sn-visible_planets_x
d2y  = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)

earth = 0

#Calculating the distance between the planets and earth
d3= planet_dist-earth

# Calculating the distance between the event and earth
d1 = event_dist-earth

# Calculating the signal time for each planetn and creating the schedule
signal_time = (d2+d3)-d1

schedule = 1987+signal_time

visible_planets['schedule'] = 1987+(signal_time)

# Filtering the planets for the specific exoplanet of interest
visible_planets = visible_planets[visible_planets['P_NAME']== 'HD 20794 d'] 

#%%
visible_planets['S_DISTANCE'] = visible_planets['S_DISTANCE']*3.26156  # changing dataset so in light years instead of parsecs

target_year = 1991
range_years = 5

# Years since the event
catalog_year = 1987  
years_since = target_year - catalog_year

# Filter the visible_stars DataFrame based on the correct range
visible_planets = visible_planets[(visible_planets['schedule']>= target_year - range_years) & 
                              (visible_planets['schedule'] <= target_year + range_years)]

# Setting up for ellipse

# Correct coordinates for targets based on the year of observation for plotting the ellipse
c_sch1 = visible_planets.loc[(visible_planets['schedule'].astype(int)) == (visible_planets['schedule'].astype(int)).max()]
c_sch  = visible_planets.loc[(visible_planets['schedule'].astype(int)) == target_year]
c_sch2 = visible_planets.loc[(visible_planets['schedule'].astype(int)) == (visible_planets['schedule'].astype(int)).min()]

# ra and dec in terms of proper motion for the stars within the target year
ra_year = c_sch['S_RA'] * years_since / 3600  # Convert arcseconds to degrees
dec_year = c_sch['S_DEC']* years_since / 3600

new_ra = c_sch['S_RA'] + ra_year
new_dec = c_sch['S_DEC'] + dec_year

new_ra = np.radians(new_ra)  
new_dec = np.radians(new_dec)

# new position of star for ellipse
new_star_xe = c_sch['S_DISTANCE'] * np.cos(new_dec) * np.cos(new_ra)
new_star_ye = c_sch['S_DISTANCE'] * np.cos(new_dec) * np.sin(new_ra)
new_star_ze = c_sch['S_DISTANCE'] * np.sin(new_dec)

# Calculate distances from the event
event_2_starx = x_sn - float(new_star_xe.values[0])
event_2_stary = y_sn - float(new_star_ye.values[0])
event_2_starz = z_sn - float(new_star_ze.values[0])
star_2_earth = float(c_sch['S_DISTANCE'].values[0])

event_2_star = np.sqrt((float(new_star_xe.values[0]) - x_sn)**2 + (float(new_star_ye.values[0]) - y_sn)**2)

# foci
sn_foci_x,sn_foci_y  = x_sn, y_sn
e_foci_x, e_foci_y = 0, 0
event_dist_2d = np.sqrt(x_sn**2+y_sn**2) # event distance in 2d 

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


fig = plt.figure()
# Plot ellipse
plt.plot(x, y, color = 'red')

plt.scatter(new_star_xe, new_star_ye, color = 'purple')
# Show focii
plt.plot(sn_foci_x, sn_foci_y, 'bo')
plt.plot(e_foci_x, e_foci_y, 'go')
plt.legend([target_year, 'stars', 'Event', 'Earth' ], fontsize = '20')
plt.xlabel('lightyears')
plt.ylabel('lightyears')
plt.axis('equal')
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

date = (visible_planets['schedule']).apply(decimal_year_to_date)

visible_planets['date'] = date

final_data = list(zip(visible_planets['S_TYPE'], visible_planets['S_RA'], visible_planets['S_DEC'], 
                      visible_planets['S_DISTANCE'], visible_planets['schedule'], visible_planets['date']))
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




