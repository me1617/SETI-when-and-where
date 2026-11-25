# -*- coding: utf-8 -*-
"""
Created on Tue May 20 10:50:58 2025

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
 #Creating event for Betelguese

Betelguese = ss.Event(ra=ss.hms2deg(5, 55, 10.30536),
                      dec=ss.deg_arcmin_arcsec2deg(7, 24, 25.4304),
                      dist=ss.pc2ly(400),
                      event_name='Betelgeuse')

# Getting Betelguese ra, dec, and distance

ra_sn = Betelguese.ra
dec_sn = Betelguese.dec
dist_sn = Betelguese.dist

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

schedule = 2030+(signal_time) # Putting all the schedule values as a whole number to simplify the process for calculating the ellipse and not choosing a specific date for the calculation
visible_stars['schedule'] = 2030+(signal_time)  # saving the schedule as the exact schedule value so when convert to proper date and time it is accurate

visible_stars['signal_time'] = signal_time

# Choosing the year we want to look at and allowing a range of 10 years on the plot
target_year = 2050
range_years = 1

# year of event and when the targets will see event
SN_year = 2030  
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
years = np.arange(target_year - range_years, target_year + range_years + 1)
SN_date = datetime.strptime('1987-02-23 05:31:00', '%Y-%m-%d %H:%M:%S')
# Current date and time
a_val = []
b_val = []
c_val = []
f_val = []
surface_val = []
for year in years:
    dt2 = datetime.strptime(f'{year}-02-23 05:31:00', '%Y-%m-%d %H:%M:%S')
    time_since_SN=ef.date_diff_in_Seconds(dt2, SN_date)
    [x_e,y_e]=ef.calc_SETI_ellipsoid(event_dist,time_since_SN,'1987a') #lyr
    L = event_dist / 2  # lyr
    T = time_since_SN * c_lyr_s  # lyr
    A = (2 * L + T) / 2  # lyr
    C= A*2
    x_e_centered = x_e-A
    E_pos=np.array([0,0])
    SN_pos=np.array([-51.4, 0])
    b= (A**2-x_e_centered**2)
    surface_area = 4*np.pi*((((A*b)**1.6075)+ ((A*C)**1.6075)+((b*C)**1.6075))/3)**(1/1.6075) # Surface area equation 
    a_val.append(A)
    b_val.append(b)
    c_val.append(C)
    surface_val.append(surface_area)
    
    fig = plt.figure()
    plt.plot(years, surface_val)
    # plt.plot(years, a_val)
    # plt.plot(years, b_val)
    # plt.plot(years, c_val)
    plt.xlabel('Years', fontsize=20)
    plt.ylabel('Surface Area', fontsize = 20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize = 20)
    plt.show()
