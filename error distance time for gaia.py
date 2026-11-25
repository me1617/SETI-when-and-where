# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 13:08:23 2025

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
votable = parse(r"C:\Users\meg16\Downloads\Data (6)\Data\GCNS_allcol.vot")  # Update with your actual file path
table = votable.get_first_table().to_table()
# Convert to Pandas DataFrame
gaia = table.to_pandas()
print(gaia.columns)


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
y_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn * np.sin(np.radians(dec_sn))


event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

d2x = x_sn-visible_planets_x
d2y  = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)

earth = 0

d3 = planet_dist-earth

d1 = event_dist-earth

signal_time = (d2+d3)-d1
signal_time = np.array(signal_time)
schedule = 1987+signal_time

visible_planets['schedule'] = 1987+np.round(signal_time)

visible_planets['signal_time'] = signal_time
visible_planets['Dist50'] = visible_planets['Dist50']*3.26156e3 # changing dataset so in light years instead of parsecs


# converting all visible stars to correct coordinates to then plot
delta_vs_ra = np.radians(visible_planets['RA_ICRS']*15) * signal_time  # Convert arcseconds to degrees
delta_vs_dec = np.radians(visible_planets['DE_ICRS']) * signal_time   # Convert arcseconds to degrees


delta_x = visible_planets['Dist50'] * signal_time * np.cos(delta_vs_dec) * np.cos(delta_vs_ra)
delta_y = visible_planets['Dist50'] * signal_time * np.cos(delta_vs_dec) * np.sin(delta_vs_ra)
delta_z = visible_planets['Dist50'] * signal_time * np.sin(delta_vs_dec)


# Recalculate the RA and Dec after applying proper motion
new_vs_ra =visible_planets['RA_ICRS'] + delta_vs_ra
new_vs_dec = visible_planets['DE_ICRS'] + delta_vs_dec

# Recalculate the Cartesian coordinates
new_vs_ra_radians = np.radians(new_vs_ra * 15)  # Convert from hours to degrees, then to radians
new_vs_dec_radians = np.radians(new_vs_dec)

new_planet_x = visible_planets['Dist50'] * np.cos(new_vs_dec_radians) * np.cos(new_vs_ra_radians)
new_planet_y = visible_planets['Dist50'] * np.cos(new_vs_dec_radians) * np.sin(new_vs_ra_radians)
new_planet_z = visible_planets['Dist50'] * np.sin(new_vs_dec_radians)

new_planet_dist = np.sqrt(new_planet_x**2+new_planet_y**2+new_planet_z**2)


new_planet_xmin = visible_planets['Dist16'] * np.cos(new_vs_dec_radians) * np.cos(new_vs_ra_radians)
new_planet_ymin = visible_planets['Dist16'] * np.cos(new_vs_dec_radians) * np.sin(new_vs_ra_radians)
new_planet_zmin = visible_planets['Dist16'] * np.sin(new_vs_dec_radians)

new_planet_dist_min = np.sqrt(new_planet_xmin**2+new_planet_ymin**2+new_planet_zmin**2)

new_planet_xmax = visible_planets['Dist84'] * np.cos(new_vs_dec_radians) * np.cos(new_vs_ra_radians)
new_planet_ymax = visible_planets['Dist84'] * np.cos(new_vs_dec_radians) * np.sin(new_vs_ra_radians)
new_planet_zmax = visible_planets['Dist84'] * np.sin(new_vs_dec_radians)

new_planet_dist_max = np.sqrt(new_planet_xmax**2+new_planet_ymax**2+new_planet_zmax**2)




#%%

planet_dist = new_planet_dist
planet_ra_new =new_vs_ra_radians
planet_dec_new = new_vs_dec_radians

dist_max_error = (new_planet_dist_max-planet_dist)/planet_dist
dist_min_error = (planet_dist-new_planet_dist_min)/planet_dist


x_sn = dist_sn * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
y_sn = dist_sn  * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
z_sn = dist_sn *np.sin(np.radians(planet_dec_new))

visible_planets_x = planet_dist * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
visible_planets_y = planet_dist * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
visible_planets_z = planet_dist * np.sin(np.radians(planet_dec_new))


max_planets_x = (planet_dist * (1 + dist_max_error)) * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
max_planets_y = (planet_dist * (1 + dist_max_error)) * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
max_planets_z = (planet_dist * (1 + dist_max_error)) * np.sin(np.radians(planet_dec_new))

min_planets_x = (planet_dist * (1 + dist_min_error)) * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
min_planets_y = (planet_dist * (1 + dist_min_error)) * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
min_planets_z = (planet_dist * (1 + dist_min_error)) * np.sin(np.radians(planet_dec_new))

max_sn = (dist_sn+(1.2e3*3.26156))
max_x_sn = max_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
max_y_sn = max_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
max_z_sn = max_sn * np.sin(np.radians(dec_sn))

max_d2x = max_x_sn - max_planets_x
max_d2y = max_y_sn - max_planets_y
max_d2z = max_z_sn - max_planets_z


max_d2 = (np.sqrt(max_d2x**2+max_d2y**2)-d2)/d2

min_sn = (dist_sn-(1.2e3*3.26156))

min_x_sn = min_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
min_y_sn = min_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
min_z_sn = min_sn * np.sin(np.radians(dec_sn))

min_d2x = min_x_sn - min_planets_x
min_d2y = min_y_sn - min_planets_y
min_d2z = min_z_sn - min_planets_z

min_d2 = (d2-np.sqrt(min_d2x**2+min_d2y**2))/d2

event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)
max_event_dist = np.sqrt(max_x_sn**2+max_y_sn**2)
min_event_dist = np.sqrt(min_x_sn**2+min_y_sn**2)


d2x = x_sn-visible_planets_x
d2y = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z



earth = 0


max_d3 = (((planet_dist * (1 + dist_max_error))+earth)-d3)/d3
min_d3 = (((planet_dist * (1 + dist_min_error))+earth)-d3)/d3


max_d1 = ((max_event_dist) - d1) / d1
min_d1 = (d1 - (min_event_dist)) / d1


min_signal_time = min_d2 + min_d3 - min_d1
max_signal_time = max_d2 + max_d3 - max_d1


# # Create subplots~
fig, ax = plt.subplots(1, 3, figsize=(12, 5))

# Min error plot
ax[0].plot(min_signal_time, dist_min_error,  label='Min Error')
ax[0].set_xlabel('Min Error Time')
ax[0].set_ylabel('Min Error Distance')
ax[0].legend()
ax[0].grid(True)

# Max error plot
ax[1].plot(max_signal_time, dist_max_error, label='Max Error', color='red')
ax[1].set_xlabel('Max Error Time')
ax[1].set_ylabel('Max Error Distance')
ax[1].legend()
ax[1].grid(True)

ax[2].plot(signal_time, planet_dist, color='red')
ax[2].legend()
ax[2].grid(True)
plt.tight_layout()
plt.show()

min_schedule = 1987+(signal_time*(1+min_signal_time))
max_schedule = 1987+(signal_time*(1+max_signal_time))

visible_planets['min_schedule'] = min_schedule
visible_planets['max_schedule'] = max_schedule

from datetime import datetime, timedelta

def decimal_year_to_date(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    days_in_year = (datetime(year + 1, 1, 1) - start_of_year).days
    min_date = start_of_year + timedelta(days=remainder * days_in_year)
    max_date = start_of_year + timedelta(days=remainder * days_in_year)
    return min_date, max_date
#%%

min_date = (visible_planets['min_schedule']).apply(decimal_year_to_date)
max_date = (visible_planets['max_schedule']).apply(decimal_year_to_date)

visible_planets['min_date'] = min_date
visible_planets['max_date'] = max_date

final_data = list(zip(visible_planets['RA_ICRS'][:30], visible_planets['DE_ICRS'][:30], visible_planets['Dist16'][:30],
                      visible_planets['Dist50'][:30],visible_planets['Dist84'][:30], visible_planets['min_schedule'][:30], visible_planets['min_date'][:30], visible_planets['max_schedule'][:30], visible_planets['max_date'][:30]))
# Column headers
headers = [ "ra", "dec","min_dist (ly)", "dist (ly)","max_dist (ly)", "min_schedule", "min_date", "max_schedule", "max_date"]


# Create a figure and axis
fig, ax = plt.subplots(figsize=(20,25))
ax.set_axis_off()  # Hide the axes

# Create the table
table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')

# Adjust font size for all cells
table.auto_set_font_size(False)
table.set_fontsize(20)  # Set your desired font size

# Optionally adjust row/column sizes
table.scale(1.2, 2)  # (col_width, row_height)
plt.tight_layout()
plt.show()

# Not a perfect linear relationship since there is fluctuations coming from the distance and time that certain planets will see the event
# However from the ellipse when coming from the time that the planets would see the signal there is some uncertainty of the exact time 
# Look at the max and min schedule next