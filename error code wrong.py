

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 14:46:41 2025

@author: meg16
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse
import setisynch.gaia as ss


# Parse the VOTable
votable = parse(r"C:\Users\meg16\Downloads\Data (2)\Data\GCNS_morecol.vot")
table = votable.get_first_table().to_table()
# Convert to Pandas DataFrame
gaia = table.to_pandas()
print(gaia.columns)


# Filter for stars visible by the Green Bank Telescope
visible_stars = gaia[(gaia['Dist50']*3.26156e3<=100000) & (gaia['DE_ICRS'] >= -46) & (gaia['DE_ICRS'] <= 90) & (gaia['RA_ICRS'] >= 0) & (gaia['RA_ICRS'] <=24)]

#Obtain the ra, dec and distance of the visible stars
ra_radians = np.radians(visible_stars['RA_ICRS'] * 15)  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_stars['DE_ICRS'])
planet_dist = visible_stars['Dist50']*3.26156e3

# Calculate the coordinates of the stars
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

# Calculating the coordinates of the event
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

visible_stars['schedule'] = 1987+np.round(signal_time)

visible_stars['signal_time'] = signal_time
visible_stars['Dist50'] = visible_stars['Dist50']*3.26156e3 # changing dataset so in light years instead of parsecs
visible_stars['Dist84'] = visible_stars['Dist84']*3.26156e3 
visible_stars['Dist16'] = visible_stars['Dist16']*3.26156e3 


target_year = 2096
range_years = 10

# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year



#%%


# Updating the ra and dec values of the stars and calculating their percentage errors and the distance percentage error
planet_ra_new = visible_stars['RA_ICRS']*15
planet_dec_new = visible_stars['DE_ICRS']

ra_error = (visible_stars['e_RA_ICRS']- visible_stars['RA_ICRS'])/visible_stars['RA_ICRS']
dec_error = (visible_stars['e_DE_ICRS'] - visible_stars['DE_ICRS'])/visible_stars['DE_ICRS']
new_ra_error = planet_ra_new*(1+ra_error) # Applying the percentage errors to the actual values
new_dec_error = planet_dec_new*(1+dec_error)

dist_max_error = (visible_stars['Dist84']-visible_stars['Dist50'])/visible_stars['Dist50']
dist_min_error = (visible_stars['Dist50']-visible_stars['Dist16'])/visible_stars['Dist50']

# Calculating the event coordinates and the stars coordinates again and then calculating their error
x_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
y_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn * np.sin(np.radians(dec_sn))

visible_planets_x = planet_dist * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
visible_planets_y = planet_dist * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
visible_planets_z = planet_dist * np.sin(np.radians(planet_dec_new))

 # When calculating cooridnates apply the percentage error to the actual distance
max_planets_x = (planet_dist * (1 + dist_max_error)) * np.cos(np.radians(new_dec_error)) * np.cos(np.radians(new_ra_error))
max_planets_y = (planet_dist * (1 + dist_max_error)) * np.cos(np.radians(new_dec_error)) * np.sin(np.radians(new_ra_error))
max_planets_z = (planet_dist * (1 + dist_max_error)) * np.sin(np.radians(new_dec_error))

min_planets_x = (planet_dist * (1 + dist_min_error)) * np.cos(np.radians(new_dec_error)) * np.cos(np.radians(new_ra_error))
min_planets_y = (planet_dist * (1 + dist_min_error)) * np.cos(np.radians(new_dec_error)) * np.sin(np.radians(new_ra_error))
min_planets_z = (planet_dist * (1 + dist_min_error)) * np.sin(np.radians(new_dec_error))

max_sn = (dist_sn+(1.2e3*3.26156))
max_x_sn = max_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
max_y_sn = max_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
max_z_sn = max_sn * np.sin(np.radians(dec_sn))

min_sn = (dist_sn-(1.2e3*3.26156))

min_x_sn = min_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
min_y_sn = min_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
min_z_sn = min_sn * np.sin(np.radians(dec_sn))

# Calculating the coordinates and distance percentage errors between the event and stars
max_d2x = max_x_sn - max_planets_x
max_d2y = max_y_sn - max_planets_y
max_d2z = max_z_sn - max_planets_z

max_d2 = (np.sqrt(max_d2x**2+max_d2y**2)-d2)/d2

min_d2x = min_x_sn - min_planets_x
min_d2y = min_y_sn - min_planets_y
min_d2z = min_z_sn - min_planets_z

min_d2 = (d2-np.sqrt(min_d2x**2+min_d2y**2))/d2

event_dist = np.sqrt(x_sn**2+y_sn**2)
max_event_dist = np.sqrt(max_x_sn**2+max_y_sn**2)
min_event_dist = np.sqrt(min_x_sn**2+min_y_sn**2)

earth = 0

# Calculating the coordinates and distance percentage errors between the stars and earth
max_d3 = (((planet_dist * (1 + dist_max_error))+earth)-d3)/d3
min_d3 = (((planet_dist * (1 + dist_min_error))+earth)-d3)/d3

# Calculating the coordinates and distance percentage errors between the event and earth
max_d1 = ((max_event_dist) - d1) / d1
min_d1 = (d1 - (min_event_dist)) / d1

# Calculating the percentage error on the signal time 
min_signal_time = min_d2 + min_d3 - min_d1
max_signal_time = max_d2 + max_d3 - max_d1

# Plotting the time and distance error relation as a subplot
fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# Min error plot
ax[0].plot(min_signal_time, dist_min_error,   label='Min Error')
ax[0].set_xlabel('Min Error Time')
ax[0].set_ylabel('Min Error Distance')
ax[0].legend()


# Max error plot
ax[1].plot(max_signal_time, dist_max_error, label='Max Error', color='red')
ax[1].set_xlabel('Max Error Time')
ax[1].set_ylabel('Max Error Distance')
ax[1].legend()



plt.tight_layout()
plt.show()

# Calculating the error on the schedule not using percentage error and saving to the file
min_schedule = 1987+(signal_time*(1+min_signal_time))
max_schedule = 1987+(signal_time*(1+max_signal_time))

visible_stars['min_schedule'] = min_schedule
visible_stars['max_schedule'] = max_schedule

from datetime import datetime, timedelta
#Creating function to convert the decimal schedule to proper date and time
def real_date(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    days_in_year = (datetime(year + 1, 1, 1) - start_of_year).days
    date = start_of_year + timedelta(days=remainder * days_in_year)
    return date

#%%

# Applying the date function to the minimum and maximum schedule for the full dataset and save them to the file
min_date = (visible_stars['min_schedule']).apply(real_date)
max_date = (visible_stars['max_schedule']).apply(real_date)

visible_stars['min_date'] = min_date
visible_stars['max_date'] = max_date

# Saving the data of interest as a zip file to create table, only containing the RA,dec, distance, schedule and date and their errors
final_data = list(zip(visible_stars['RA_ICRS'], visible_stars['DE_ICRS'], visible_stars['Dist16'],
                      visible_stars['Dist50'],visible_stars['Dist84'], visible_stars['min_schedule'],
                      visible_stars['min_date'], visible_stars['max_schedule'], visible_stars['max_date']))
# Column headers
headers = [ "ra", "dec","min_dist (ly)", "dist (ly)","max_dist (ly)", "min_schedule", "min_date", "max_schedule", "max_date"]



fig, ax = plt.subplots(figsize=(20,25))
ax.set_axis_off()  # Hide the axes

# Create a table
table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')


# Adjust font size and column width
table.auto_set_font_size(False)
table.set_fontsize(10)
table.auto_set_column_width([i for i in range(len(headers))])
