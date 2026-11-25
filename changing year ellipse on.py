# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:15:39 2025

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
target_year = 2026
range_years = 1
for year in range(2025, 2030 + 1):
    # year of event and when the targets will see event
    SN_year = 1987  
    years_since = year - SN_year
    
    # Filtering the catalogue based on the target year
    visible_stars_filtered = visible_stars[(visible_stars['schedule'] >= year - range_years) & 
                                  (visible_stars['schedule'] <= year + range_years)]
    
    print(len(visible_stars))
    print(len(visible_stars_filtered))
    
    #print(visible_stars_filtered['schedule'])
    
    
    c_m_s = 3e8
    c_lyr_s = c_m_s/9.461e+15
    
    SN_date = datetime.strptime('1987-02-23 05:31:00', '%Y-%m-%d %H:%M:%S')
    # Current date and time
    #dt2 = datetime.now()
    dt2 = datetime.strptime(f'{year}-02-23 05:31:00', '%Y-%m-%d %H:%M:%S')
    time_since_SN=ef.date_diff_in_Seconds(dt2, SN_date)
    #time_since_SN6=ef.date_diff_in_Seconds(dt6, SN_date)
    [x_e,y_e]=ef.calc_SETI_ellipsoid(event_dist,time_since_SN,'1987a') #lyr
    L = event_dist / 2  # lyr
    T = time_since_SN * c_lyr_s  # lyr
    A = (2 * L + T) / 2  # lyr
    
    x_e_centered = x_e-A
    E_pos=np.array([0,0])
    SN_pos=np.array([-51.4, 0])


    # [x_e,y_e]=ef.calc_SETI_ellipsoid(z_y_dist,time_since_SN,'1987a') #lyr
    fig = plt.figure()    # z_max = np.max(visible_stars_filtered.zcoord50)
    z_max = np.max(visible_stars_filtered.zcoord50)
    z_min = np.min(visible_stars_filtered.zcoord50)
    y_max = np.max(visible_stars_filtered.ycoord50)
    y_min = np.min(visible_stars_filtered.ycoord50)
    z_mean = np.mean(visible_stars_filtered.zcoord50)
    y_mean = np.mean(visible_stars_filtered.ycoord50)
    z_y_dist_max  = np.sqrt(z_max**2 + y_max**2)
    z_y_dist_min = np.sqrt(z_min**2 + y_min**2)
    z_y_dist = z_y_dist_max-z_y_dist_min
    [x_e,y_e]=ef.calc_SETI_ellipsoid(event_dist,time_since_SN,'1987a')
    # Create a grid over the extremes of the z and y coordinates
    z = np.linspace(z_min, z_max, 100)
    y = np.linspace(y_min, y_max, 100)
    Z, Y = np.meshgrid(z, y)
    d1 = np.sqrt((Z - E_pos[0])**2 + (Y - E_pos[1])**2)
    d2 = np.sqrt((Z - SN_pos[0])**2 + (Y - SN_pos[1])**2)
    sum_distances = d1 + d2
   # levels = np.linspace(z_max,  + 1, 5)


    #plt.scatter(visible_stars_filtered.zcoord50, visible_stars_filtered.ycoord50,c=visible_stars_filtered.xcoord50,s=1)
    plt.scatter(visible_stars_filtered.zcoord50, visible_stars_filtered.ycoord50,s=1,c=visible_stars_filtered.xcoord50,  # color by z-coordinate
        cmap='viridis'  )
    # contours = plt.contour(Z, Y, sum_distances, levels=z, colors='orange', linewidths=1)
    # plt.clabel(contours, inline=True, fontsize=10, fmt='%.1f')
    plt.plot(z_max, y_max)
    plt.plot(z_min, y_min)
    #contours = plt.contour(Z, Y, scalar_field, levels=10, colors='red', linewidths=0.5, label=f'Year {year}')
    #plt.clabel(contours, inline=True, fontsize=10, fmt='%.0f')
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
fig = plt.figure(figsize=(10, 8))

for year in range(2025, 2030 + 1):
    # Year of SN event
    SN_year = 1987  
    years_since = year - SN_year
    
    # Filter stars visible in ±1 year around the current year
    visible_stars_filtered = visible_stars[(visible_stars['schedule'] >= year - 1) & 
                                           (visible_stars['schedule'] <= year + 1)]

    print(f"Year {year}: {len(visible_stars_filtered)} stars")

    # Plot stars for this year
    plt.scatter(
        visible_stars_filtered.zcoord50, 
        visible_stars_filtered.ycoord50, s=2,
        label=f'{year}',  # Label for the legend
        alpha=0.6
    )
    plt.plot(x_e, y_e)
# Optional: Show Earth and SN positions (scaled)
plt.plot(0, 0, 'go', label='Earth')

# Plot settings
plt.xlabel('Z (lightyears)', fontsize=16)
plt.ylabel('Y (lightyears)', fontsize=16)
plt.xlim([-200, 200])
plt.ylim([-200, 200])
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(title='Observation Year', fontsize=12)
plt.title('Overlapping Star Observations by Year', fontsize=18)
plt.grid(True)
plt.tight_layout()
plt.show()
#%%

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

for year in range(2025, 2030 + 1):
    visible_stars_filtered = visible_stars[
        (visible_stars['schedule'] >= year - 1) &
        (visible_stars['schedule'] <= year + 1)
    ]

    z = visible_stars_filtered['zcoord50'].dropna()
    y = visible_stars_filtered['ycoord50'].dropna()

    if len(z) == 0 or len(y) == 0:
        continue  # Skip if data is empty

    # Compute extremes
    z_min, z_max = z.min(), z.max()
    y_min, y_max = y.min(), y.max()

    # Compute center and widths
    width = z_max - z_min
    height = y_max - y_min
    center = ((z_max + z_min) / 2, (y_max + y_min) / 2)

    # Draw ellipse
    ellipse = Ellipse(xy=center, width=width, height=height,
                      edgecolor='C0', facecolor='none', lw=2, label=f'{year}')
    plt.gca().add_patch(ellipse)

# Optional scatter to show stars
plt.scatter(visible_stars['zcoord50'], visible_stars['ycoord50'], s=1, alpha=0.1)

plt.xlabel('lightyears')
plt.ylabel('lightyears')
plt.legend()
plt.xlim(-200, 200)
plt.ylim(-200, 200)
plt.gca().set_aspect('equal')
plt.show()

    # z_max = np.max(visible_stars_filtered.zcoord50)
    # z_min = np.min(visible_stars_filtered.zcoord50)
    # y_max = np.max(visible_stars_filtered.ycoord50)
    # y_min = np.min(visible_stars_filtered.ycoord50)
    # z_mean = np.mean(visible_stars_filtered.zcoord50)
    # y_mean = np.mean(visible_stars_filtered.ycoord50)
    # z_y_dist_max  = np.sqrt(z_max**2 + y_max**2)
    # z_y_dist_min = np.sqrt(z_min**2 + y_min**2)
    # z_y_dist = z_y_dist_max-z_y_dist_min
    # # Create a grid over the extremes of the z and y coordinates
    # z = np.linspace(z_min, z_max, 300)
    # y = np.linspace(y_min, y_max, 300)
    # Z, Y = np.meshgrid(z, y)
    
    # # Define a scalar field — example: distance from (z_mean, y_mean)
    # scalar_field = np.sqrt((Z - z_mean)**2 + (Y - y_mean)**2)
    