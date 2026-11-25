# -*- coding: utf-8 -*-
"""
Created on Fri May 16 19:08:44 2025

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
range_years = 0.1

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
z_mean = np.mean(visible_stars_filtered.zcoord50)
y_mean  = np.mean(visible_stars_filtered.ycoord50)
z_max = np.max(visible_stars_filtered.zcoord50)
z_min = np.min(visible_stars_filtered.zcoord50)
y_max = np.max(visible_stars_filtered.ycoord50)
y_min = np.min(visible_stars_filtered.ycoord50)

# Foci

#event_dist_2d = np.sqrt(z_sn**2+y_sn**2) # event distance in 2d 
event_2_star_x = 0 - visible_stars_filtered.xcoord50
event_2_star_y = SN_pos[1]*3261.56 - visible_stars_filtered.ycoord50
event_2_star_z = SN_pos[0]*3261.56 - visible_stars_filtered.zcoord50

event_2_star = np.sqrt(event_2_star_x**2+event_2_star_y**2+ event_2_star_z**2)
star_2_earth = visible_stars_filtered.Dist50*3261.56
sn_foci_x,sn_foci_y  = z_sn, -y_sn
# Ellipse parameters
x0 = (SN_pos[0]*3261.56 + 0) / 2                      # Center x-value
y0 = (SN_pos[1]*3261.56+ 0) / 2                       # Center y-value
f = np.sqrt((SN_pos[0]*3261.56 - x0)**2 + (SN_pos[1]*3261.56 - y0)**2) # Distance from center to focus

c = event_2_star+star_2_earth

a = A /2                         # Semimajor axis
b = np.sqrt(A**2 - f**2)                 # Semiminor axis
phi = np.arctan2((0 - SN_pos[1]*3261.56), (0 - SN_pos[0]*3261.56))   # Angle between major y-axis and x-axis

resolution = 10000
t = np.linspace(0, 2*np.pi, resolution)

# Coordinates for the ellipse
x = x0 + A * np.cos(t) * np.cos(phi) - b * np.sin(t) * np.sin(phi)
y = y0 + A * np.cos(t) * np.sin(phi) + b * np.sin(t) * np.cos(phi)

fig = plt.figure()
# Plot ellipse
plt.plot(x, y, color = 'red')
# Plot foci
plt.plot(E_pos* 3261.56, 'go')
plt.plot(SN_pos[0]*3261.56, -SN_pos[0]*3261.56,  'ro')
plt.legend([2025, 'Planets', 'Event', 'Earth' ], fontsize = '20')
plt.xlabel('lightyears')
plt.ylabel('lightyears')
plt.axis('equal')


#plt.scatter(visible_stars_filtered.zcoord50, visible_stars_filtered.ycoord50,c=visible_stars_filtered.xcoord50,s=1)
plt.scatter(visible_stars_filtered.zcoord50, visible_stars_filtered.ycoord50,s=1,c=visible_stars_filtered.xcoord50,  # color by z-coordinate
    cmap='viridis'  )
#plt.scatter(z_max, y_max)
#plt.scatter(z_min, y_min)
# Show focii


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
z_new =[]
y_new = []
for year in range(2025, 2030+1):
    SN_date = datetime.strptime('1987-02-23 05:31:00', '%Y-%m-%d %H:%M:%S')
    # Current date and time
    dt2 = datetime.strptime(f'{year}-02-23 05:31:00', '%Y-%m-%d %H:%M:%S')
    time_since_SN=ef.date_diff_in_Seconds(dt2, SN_date)
    [x_e,y_e]=ef.calc_SETI_ellipsoid(event_dist,time_since_SN,'1987a') #lyr
    L = event_dist / 2  # lyr
    T = time_since_SN * c_lyr_s  # lyr
    A = (2 * L + T) / 2  # lyr
    # --- Define Foci in z-y plane ---
    z_earth, y_earth = 0, 0
    z_sn = dist_sn * np.sin(np.radians(dec_sn))  # z already computed earlier
    y_sn = dist_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))  # consistent with y from earlier
    
    # Foci positions
    f1 = np.array([z_earth, y_earth])
    f2 = np.array([z_sn, y_sn])
    
    # Center of the ellipse (midpoint between foci)
    center = (f1 + f2) / 2
    z0, y0 = center
    
    # Distance from center to each focus
    f = np.linalg.norm(f1 - center)
    
    # Semimajor and semiminor axes
    a = A   # total light distance from event to Earth
    b = np.sqrt(a**2 - f**2)
    
    # Angle of major axis (rotation in z-y plane)
    phi = np.arctan2(y_sn - y_earth, z_sn - z_earth)
    
    # Generate ellipse in z-y
    t = np.linspace(0, 2 * np.pi, 1000)
    z = z0 + a * np.cos(t) * np.cos(phi) - b * np.sin(t) * np.sin(phi)
    y = y0 + a * np.cos(t) * np.sin(phi) + b * np.sin(t) * np.cos(phi)
    z_new.append(z)
    y_new.append(y)
    
    # --- Plot ---
 # Calculating coordinates for the event
ra_radians = np.radians((visible_stars_filtered['RA_ICRS']+(visible_stars_filtered['RA_ICRS']*time_since_SN)/3600))  # degrees to radians
dec_radians = np.radians((visible_stars_filtered['DE_ICRS']+(visible_stars_filtered['DE_ICRS']*time_since_SN)/3600))


star_dist = visible_stars_filtered['Dist50']*3.26156e3  # extracting the distance values and converting units from kpc to light years

#Calculating the x,y,z position of the stars
visible_stars_x = star_dist * np.cos(dec_radians) * np.cos(ra_radians)
visible_stars_y = star_dist  * np.cos(dec_radians) * np.sin(ra_radians)
visible_stars_z = star_dist *np.sin(dec_radians)
fig = plt.figure(figsize=(12, 8))
plt.scatter(visible_stars_z, visible_stars_y,
            s=1, c=visible_stars_x, cmap='viridis', label='Planets')

plt.plot(z_new, y_new, 'r-', label='2025 Ellipsoid')
plt.plot([z_earth], [y_earth], 'go', label='Earth')
plt.plot([z_sn], [y_sn], 'ro', label='SN1987A')

plt.xlabel('lightyears (z)', fontsize=16)
plt.ylabel('lightyears (y)', fontsize=16)
plt.legend(fontsize=14)
plt.colorbar(label='xcoord50')
plt.axis('equal')
#plt.xlim([-200, 200])
#plt.ylim([-200, 200])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title('SETI Ellipsoid (z vs y plane)', fontsize=18)
plt.show()
#%%

from sklearn.decomposition import PCA
from matplotlib.patches import Ellipse

# Extract the 2D coordinates
coords = np.vstack((visible_stars_filtered.zcoord50, visible_stars_filtered.ycoord50)).T

# Perform PCA
pca = PCA(n_components=2)
pca.fit(coords)

# Get ellipse parameters
center = pca.mean_
width, height = 2 * np.sqrt(pca.explained_variance_)  # 1 std dev
angle = (np.arctan2(pca.components_[0,1], pca.components_[0,0]))

# Plot
fig, ax = plt.subplots(figsize=(12, 8))
sc = ax.scatter(visible_stars_filtered.zcoord50, visible_stars_filtered.ycoord50,
                s=1, c=visible_stars_filtered.xcoord50, cmap='viridis', label='Planets')

# Add fitted ellipse
ellipse = Ellipse(xy=center, width=width*2, height=height*2, angle=angle,
                  edgecolor='red', fc='None', lw=2, label='Fitted Ellipse')
ax.add_patch(ellipse)

# Decorations
ax.set_xlabel('lightyears (z)', fontsize=16)
ax.set_ylabel('lightyears (y)', fontsize=16)
ax.legend(fontsize=14)
plt.colorbar(sc, label='xcoord50')
ax.set_xlim([-200, 200])
ax.set_ylim([-200, 200])
ax.set_aspect('equal')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title('Fitted Ellipse to Star Distribution', fontsize=18)
plt.show()
