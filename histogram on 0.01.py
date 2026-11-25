# -*- coding: utf-8 -*-
"""
Created on Wed May 14 10:57:08 2025

@author: meg16
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse
import setisynch.gaia as ss
import pandas as pd


star_data = pd.read_csv(r"C:\Users\meg16\Downloads\Data (6)\Data\on_el_th0.01lyr.csv")

# Filter for stars visible by the Green Bank Telescop~e
visible_stars = star_data[(star_data['Dist50']*3.26156e3<=100000) & (star_data['DE_ICRS'] >= -46) & (star_data['DE_ICRS'] <= 90) & (star_data['RA_ICRS'] >= -360) & (star_data['RA_ICRS'] <=360)]

ra_radians = np.radians(visible_stars['RA_ICRS'])  # degrees to radians
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
event_dist = np.sqrt(x_sn**2+y_sn**2)

# Calculating the distance between the event and the stars
d2x = x_sn-visible_stars_x
d2y  = y_sn-visible_stars_y
d2z = z_sn-visible_stars_z

d2 = np.sqrt(d2x**2+d2y**2)

earth = 0

# Calculating the distance between the stars and the earth
d3 = star_dist-earth

# Calculating the distande between the event and earth
d1 = event_dist-earth

# Calculating the signal time and the schedule and saving it to the dataset
signal_time = (d2+d3)-d1
schedule = 1987+(signal_time) # Putting all the schedule values as a whole number to simplify the process for calculating the ellipse and not choosing a specific date for the calculation
visible_stars['schedule'] = 1987+(signal_time)  # saving the schedule as the exact schedule value so when convert to proper date and time it is accurate

visible_stars['signal_time'] = signal_time
visible_stars['Dist50'] = visible_stars['Dist50']*3.26156e3 # changing dataset so in light years instead of kpc
visible_stars['Dist84'] = visible_stars['Dist84']*3.26156e3 
visible_stars['Dist16'] = visible_stars['Dist16']*3.26156e3 

# Choosing the year we want to look at and allowing a range of 10 years on the plot
target_year = 2025
range_years = 10

# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year


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
star_dist_min = (visible_stars['Dist50'])*3.26156e3  # extracting the distance values and converting units from kpc to light years
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
import statistics
import pandas as pd
#Calculating the mean, mode and the median of the dataset
min_signal_years = schedule+(signal_time*(1+signal_time_min))
max_signal_years = schedule-(signal_time*(1+signal_time_max))

mean_time = np.mean(schedule)
mode_stars = statistics.mode(schedule)
median_time = np.median(schedule)
no_stars = len(visible_stars)
# Calculate the min and max years
# minmin_year = (np.min(min_signal_years))
# maxmin_year = (np.max(min_signal_years))
# minmax_year = (np.min(max_signal_years))
# maxmax_year = (np.max(max_signal_years))
# # Create year bins
# bins_min = np.arange(minmin_year, maxmin_year + 2)  # +2 to account for the last year 
# bins_max = np.arange(minmax_year, maxmax_year + 2) 
# Plot histogram
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
#plt.hist(schedule, bins=bins, color='orange')
ax[0].hist(min_signal_years, 700, color = 'orange')
ax[0].set_ylabel('No. of Stars', fontsize = '20')
ax[0].set_xlabel('Minimum Error on Year', fontsize = '20')
ax[0].tick_params(axis='x', labelsize=20)
ax[0].tick_params(axis='y', labelsize=20)
#plt.hist(schedule, bins=bins, color='orange')
ax[1].hist(max_signal_years, 700, color = 'blue')
# plt.axvline(mean_time, color = 'red')
#plt.axvline(mode_stars, color = 'blue')
#plt.axvline(median_time, color = 'green')
#plt.legend(['mean Year', 'Median Year' ])
ax[1].set_ylabel('No. of Stars', fontsize = '20')
ax[1].set_xlabel('Maximum Error on Year', fontsize = '20')
ax[1].tick_params(axis='x', labelsize=20)
ax[1].tick_params(axis='y', labelsize=20)

# Get the histogram counts and bin edges
counts_min, bin_edges_min = np.histogram(min_signal_years, bins=700)
counts_max, bin_edges_max = np.histogram(max_signal_years, bins=700)
