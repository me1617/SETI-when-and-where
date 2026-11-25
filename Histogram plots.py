# -*- coding: utf-8 -*-
"""
Created on Sun May 18 19:07:35 2025

@author: meg16
"""

import setisynch.gaia as ss
import sys
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


votable = parse(r"C:\Users\meg16\Downloads\Data (6)\Data\GCNS_allcol.vot") # opening GCNS as vot 
table = votable.get_first_table().to_table() # converting vot to table


gaia = table.to_pandas() # converting table to panda format 
print(gaia.columns)
#%%

# Filter for stars visible by the Green Bank Telescope
visible_stars = gaia[(gaia['Dist50']*3.26156e3<=100000) & (gaia['DE_ICRS'] >= -46) & (gaia['DE_ICRS'] <= 90) & (gaia['RA_ICRS'] >= -360) & (gaia['RA_ICRS'] <=360)].copy()

ra_radians = np.radians(visible_stars['RA_ICRS'])
  # degrees to radians
dec_radians = np.radians(visible_stars['DE_ICRS'])


star_dist = visible_stars['Dist50']*3.26156e3  # extracting the distance values and converting units from kpc to light years

#Calculating the x,y,z position of the stars
visible_stars_x = star_dist * np.cos(dec_radians) * np.cos(ra_radians)
visible_stars_y = star_dist  * np.cos(dec_radians) * np.sin(ra_radians)
visible_stars_z = star_dist *np.sin(dec_radians)

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
d2x = x_sn-visible_stars_x
d2y  = y_sn-visible_stars_y
d2z = z_sn-visible_stars_z

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)
#print(d2[0])
earth = 0

# Calculating the distance between the stars and the earth
d3 = star_dist-earth

# Calculating the distande between the event and earth

d1 = event_dist
#print(event_dist)

# Calculating the signal time and the schedule and saving it to the dataset
signal_time = (d2+d3)-d1 # EC added as I think it should be in 3D
print(signal_time[0])

schedule = 1987+(signal_time) # Putting all the schedule values as a whole number to simplify the process for calculating the ellipse and not choosing a specific date for the calculation
visible_stars['schedule'] = 1987+(signal_time)  # saving the schedule as the exact schedule value so when convert to proper date and time it is accurate

visible_stars['signal_time'] = signal_time

# Choosing the year we want to look at and allowing a range of 10 years on the plot
target_year = 2025
range_years = 1

# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year

#%%
# Filtering the catalogue based on the target year
visible_stars_filtered = visible_stars[(visible_stars['schedule'] >= target_year - range_years) & 
                              (visible_stars['schedule'] <= target_year + range_years)]

print(len(visible_stars))
print(len(visible_stars_filtered))

#print(visible_stars_filtered['schedule'])

#%%

c_m_s = 3e8
c_lyr_s = c_m_s/9.461e+15

SN_date = datetime.strptime('1987-02-23 05:31:00', '%Y-%m-%d %H:%M:%S')
# Current date and time
dt2 = datetime.now()
time_since_SN=ef.date_diff_in_Seconds(dt2, SN_date)
[x_e,y_e]=ef.calc_SETI_ellipsoid(event_dist,time_since_SN,'1987a') #lyr
L = event_dist / 2  # lyr
T = time_since_SN * c_lyr_s  # lyr
A = (2 * L + T) / 2  # lyr

x_e_centered = x_e-A
E_pos=np.array([0,0])
SN_pos=np.array([-51.4, 0])

fig = plt.figure()

#plt.scatter(visible_stars_filtered.zcoord50, visible_stars_filtered.ycoord50,c=visible_stars_filtered.xcoord50,s=1)
plt.scatter(visible_stars_filtered.zcoord50, visible_stars_filtered.ycoord50,s=1,c=visible_stars_filtered.xcoord50,  # color by z-coordinate
    cmap='viridis'  )

# Show focii
plt.plot(E_pos* 3261.56, 'go')
plt.xlabel('lightyears', fontsize = '20')
plt.ylabel('lightyears', fontsize = '20')
#plt.axis('equal')
cbar = plt.colorbar(label='xcoord50')
cbar.ax.tick_params(labelsize=20)      # Adjusting the fontsize of the colourbar
cbar.set_label('xcoord50', fontsize=20)

plt.ylim([-200, 200])
plt.xlim([-200, 200])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()

#%%

#Calculation for the errors in ra, dec and distance
star_ra_new = np.radians(visible_stars['RA_ICRS'])
star_dec_new = np.radians(visible_stars['DE_ICRS'])


#Calculating the percentage errors for ra, dec and distance EC: x 100 if percentage error
ra_error = ((visible_stars['e_RA_ICRS'])/(3.6e6))/(visible_stars['RA_ICRS']) 
dec_error = ((visible_stars['e_DE_ICRS'])/(3.6e6))/visible_stars['DE_ICRS']

dist_max_error = visible_stars['Dist84']/visible_stars['Dist50']
dist_min_error = visible_stars['Dist16']/visible_stars['Dist50']

# Applying the percentage error to the actual values of ra and dec # EC: Isn't this multiplying by a radian value?
new_ra_error = star_ra_new*(1+ra_error)
new_dec_error = star_dec_new*(1+dec_error)


# x_sn = dist_sn * np.cos(np.radians(new_dec_error)) * np.cos(np.radians(new_ra_error))
# y_sn = dist_sn  * np.cos(np.radians(new_dec_error)) * np.sin(np.radians(new_ra_error))
# z_sn = dist_sn *np.sin(np.radians(new_dec_error))

# visible_stars_x = star_dist * np.cos(np.radians(new_dec_error)) * np.cos(np.radians(new_ra_error))
# visible_stars_y = star_dist * np.cos(np.radians(new_dec_error)) * np.sin(np.radians(new_ra_error))
# visible_stars_z = star_dist * np.sin(np.radians(new_dec_error))


def error_coordinates(dist_error ,ra, dec, ra_e, dec_e):
    cos_r_error = -np.sin(np.radians(ra))*np.radians(ra_e) # cosine error propagating equation
    cos_r_actual = np.cos(np.radians(ra))
    cos_d_error = -np.sin(np.radians(dec))*np.radians(dec_e) # cosine error propagating equation
    cos_d_actual = np.cos(np.radians(dec))
    sin_r_error = np.cos(np.radians(ra))*np.radians(ra_e) # sine error propagating equation
    sin_r_actual = np.sin(np.radians(ra))
    sin_d_error = np.cos(np.radians(dec))*np.radians(dec_e) # sine error propagating equation
    sin_d_actual = np.sin(np.radians(dec))
    x_error = np.sqrt((dist_error)**2 + (cos_d_error/cos_d_actual)**2 + (cos_r_error/cos_r_actual)**2)
    y_error = np.sqrt((dist_error)**2 + (cos_d_error/cos_d_actual)**2 + (sin_r_error/sin_r_actual)**2)
    z_error = np.sqrt((dist_error)**2 + (sin_d_error/sin_d_actual)**2)
    return x_error, y_error, z_error
    

def distance(x_1, y_1, z_1, x_2, y_2, z_2):
    return np.sqrt((x_1-x_2)**2+(y_1-y_2)**2+(z_1-z_2)**2)

#Whole number
# dist_max = (star_dist * (1 + dist_max_error))
# dist_min = (star_dist * (1 + dist_min_error))

#Decimal 
max_stars_x,max_stars_y,max_stars_z = error_coordinates(dist_max_error,star_ra_new, star_dec_new, ra_error, dec_error) 
min_stars_x, min_stars_y ,min_stars_z  = error_coordinates(dist_min_error,star_ra_new, star_dec_new, ra_error, dec_error)

#Whole number
d2 = distance(x_sn, y_sn,z_sn, visible_stars_x, visible_stars_y, visible_stars_z)

# Calculating the error maximum distance between the event and the stars
max_d2 = (distance(x_sn, y_sn, z_sn, max_stars_x*visible_stars_x, max_stars_y*visible_stars_y, max_stars_z*visible_stars_z))

# Calculating the minimum percentage error distance between the event and the stars
min_d2 = (distance(x_sn, y_sn, z_sn,  min_stars_x*visible_stars_x, min_stars_y*visible_stars_y, min_stars_z*visible_stars_z))

event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)
max_event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)
min_event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

earth = 0
d3 = star_dist-earth
d1 = event_dist-earth

# Calculating the percentage error distance between the stars and earth
max_d3 = (((visible_stars['Dist84']*3.26156e3)-earth))
min_d3 = (((visible_stars['Dist16']*3.26156e3)-earth))

# Calculating the percentage error distance between the event and earth
max_d1 = ((max_event_dist)) 
min_d1 =  (min_event_dist) 

signal_time = (d2+d3)-d1
# # Calculating the signal_time error using prop
# min_signal_time =np.sqrt(min_d2**2+min_d3**2+min_d1**2)#/signal_time
# max_signal_time =np.sqrt(max_d2**2+max_d3**2+max_d1**2)#/signal_time


#%%
# Error Calculation. Rerun for the error limits to get time limits.

star_dist_max = (visible_stars['Dist84'])*3.26156e3  # extracting the distance values and converting units from kpc to light years
star_dist_min = (visible_stars['Dist16'])*3.26156e3  # extracting the distance values and converting units from kpc to light years
# Calculating the distance between the stars and the earth
d3_max = star_dist_max-earth
d3_min = star_dist_min-earth

# Calculating the signal time and the schedule and saving it to the dataset
signal_time_max = (max_d2+max_d3)-max_d1 # EC added as I think it should be in 3D
signal_time_min = (min_d2+min_d3)-min_d1 # EC added as I think it should be in 3D


# Creating the minimum and maximum schedule
# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year

# Filtering the catalogue based on the target year
# visible_stars_filtered = visible_stars[(visible_stars['schedule'] >= target_year - range_years) & 
#                               (visible_stars['schedule'] <= target_year + range_years)]

print(visible_stars['schedule'])

# So you could filter within the whole error window of each star... and you could perhaps plot the error as a colourbar on the ellipse plot.

#%%
# Ensuring the signal errors are applied to the correct filtered data
filtered_error_max = signal_time_max[visible_stars.index]
filtered_error_min = signal_time_min[visible_stars.index]
# Calculating the error on the schedule 
min_schedule = 1987+filtered_error_min
max_schedule = 1987+filtered_error_max
# Calculating the difference in signal error value to be the time window value 
difference = filtered_error_max-filtered_error_min
now_difference_max = filtered_error_max-signal_time[visible_stars.index]
visible_stars['difference'] = difference # Saving the time window values to the file
filter_time_window = visible_stars_filtered[visible_stars['difference']<2.8] # Filtering the candidates based off the mean time window value
# Saving the error values for the schedule to the file
filter_time_window['min_schedule']  = min_schedule
filter_time_window['max_schedule']  = max_schedule
# Plotting the new filtered candidates in the z-y plane with the time window as the colourbar
fig = plt.figure()

#plt.scatter(visible_stars_filtered.zcoord50, visible_stars_filtered.ycoord50,c=visible_stars_filtered.xcoord50,s=1)
plt.scatter(filter_time_window.zcoord50, filter_time_window.ycoord50,s=3,c=(filtered_error_max[filter_time_window.index]-
                                                                            filtered_error_min[filter_time_window.index]),  # color by z-coordinate
    cmap='viridis'  )

# Show focii
plt.plot(E_pos* 3261.56, 'go')

#plt.axis('equal')
cbars = plt.colorbar(label='Time Window (Years)')
cbars.ax.tick_params(labelsize=20)      # Adjusting the fontsize of the colourbar
cbars.set_label('Time Window (Years)', fontsize=20)
plt.xlabel('Z (Lightyears)', fontsize = 20)
plt.ylabel('Y (Lightyears)', fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.ylim([-200,200])
plt.xlim([-200,200])
plt.show()
print(max_schedule)
print(np.mean(max_schedule))
print(min_schedule)
print(np.mean(min_schedule))
#%%

from datetime import datetime, timedelta
#Creating function to convert the decimal schedule to proper date and time
def real_date(decimal_year):
    year = int(decimal_year)
    remainder = decimal_year - year
    start_of_year = datetime(year, 1, 1)
    days_in_year = (datetime(year + 1, 1, 1) - start_of_year).days
    date = start_of_year + timedelta(days=remainder * days_in_year)
    return date

# Extracting the first 30 rows for the c_sch_min calculation and applying the time function to the original and minimum schedule
min_date = (filter_time_window['min_schedule']).apply(real_date)
max_date = (filter_time_window['max_schedule']).apply(real_date)
date = (filter_time_window['schedule']).apply(real_date)

# Then saving them in the file 
filter_time_window['date'] = date
filter_time_window['min_date'] = min_date
filter_time_window['max_date'] = max_date
# Placing the information of interest in a zip file to then plot as a table. Rounded to 3 d.p. for the first 10 candidates
final_data = list(zip( filter_time_window[:10]['GaiaEDR3'], round(filter_time_window[:10]['RA_ICRS'], 3), round(filter_time_window[:10]['DE_ICRS'], 3),
                      round(filter_time_window[:10]['Dist50']*3261.56, 3),round(filter_time_window[:10]['date'], 3), round(filter_time_window[:10]['difference']*12, 3)))

# Creating the table headers and the table
headers = [ "GaiaEDR3", "RA (Degrees)", "Dec (Degrees)", "Distance (Lyrs)", "Date (Y/M/D, H:M:S)", "Time Window (Months)"]

fig, ax = plt.subplots(figsize=(30, 35)) 
ax.set_axis_off()

table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')

# Adjusting the size of the table cells
table.auto_set_font_size(False)
table.set_fontsize(15)  

table.scale(0.8, 2)
plt.tight_layout()
plt.show()
#%%
# Histogram of the full schedule of all stars seen with GBT
fig = plt.figure()
#plt.hist(schedule, bins=bins, color='orange')
plt.hist(schedule, 700, color = 'orange')
plt.ylabel('No. of Stars', fontsize = '20')
plt.xlabel('Years', fontsize = '20')
plt.xticks( fontsize=20)
plt.yticks(fontsize=20)
# Get the histogram counts and bin edges
counts_min, bin_edges_min = np.histogram(schedule, bins=700)
#counts_max, bin_edges_max = np.histogram(max_signal_years, bins=700)