# -*- coding: utf-8 -*-
"""
Created on Tue May 20 18:37:22 2025

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

#Creating event for Betelguese

Betelguese = ss.Event(ra=ss.hms2deg(5, 55, 10.30536),
                      dec=ss.deg_arcmin_arcsec2deg(7, 24, 25.4304),
                      dist=ss.pc2ly(400),
                      event_name='Betelgeuse')

# Getting Betelguese ra, dec, and distance

ra_sn = Betelguese.ra
dec_sn = Betelguese.dec
dist_sn = Betelguese.dist

x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn *np.sin(np.radians(dec_sn))


event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

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

schedule = 2300+signal_time

visible_planets['schedule'] = 2300+(signal_time)

# Filtering the planets for the specific exoplanet of interest
# visible_HD2 = visible_planets['P_NAME'].str.startswith('HD 219134')
# visible_planets = visible_planets[visible_planets['P_NAME'].isin([ 'TOI-1693 b', 'HD 33564 b'])]
# visible_planets = pd.merge(visible_HD2, visible_planets)
# Mask for all planets in the HD 219134 system
mask_hd = visible_planets['P_NAME'].str.startswith('HD 219134')

# Mask for the other specific planets
mask_specific = visible_planets['P_NAME'].isin(['TOI-1693 b', 'HD 33564 b'])

# Combine both masks
visible_planets = visible_planets[mask_hd | mask_specific]


#%%
visible_planets['S_DISTANCE'] = visible_planets['S_DISTANCE']*3.26156  # changing dataset so in light years instead of parsecs


#converting the decimal year to specific date

from datetime import datetime, timedelta

def decimal_year_to_date(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    days_in_year = (datetime(year + 1, 1, 1) - start_of_year).days
    date = start_of_year + timedelta(days=remainder * days_in_year)
    return date.date()


date = (visible_planets['schedule']).apply(decimal_year_to_date)

visible_planets['date'] = date

final_data = list(zip(visible_planets['P_NAME'], visible_planets['S_RA'], visible_planets['S_DEC'], 
                      visible_planets['S_DISTANCE'], visible_planets['schedule'], visible_planets['date']))
# Column headers
headers = ["Name", "RA (Degrees)", "Dec (Degrees)", "Distance (ly)", "Schedule", "Date (Y/M/D)"]

# Creating the table
fig, ax = plt.subplots(figsize=(5, 3))
ax.set_axis_off() 

table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')

# Set font size
table.auto_set_font_size(False)  
table.set_fontsize(20)  # Change the size as needed



# Show the figure
plt.show()




