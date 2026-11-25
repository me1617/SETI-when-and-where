# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 18:15:02 2025

@author: meg16
"""

from astropy.io.votable import parse
import pandas as pd
#matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import setisynch.gaia as ss
import setisynch.hygdata_v3 as hs

# Parse the VOTable
votable = parse(r"C:\Users\meg16\Downloads\Data (1)\Data\GCNS.vot")  # Update with your actual file path
table = votable.get_first_table().to_table()
 
# Convert to Pandas DataFrame
stars_gaia = table.to_pandas()
 
ra_hyg = np.asanyarray(stars_gaia['RA_ICRS']) # in degrees
dec_hyg = np.asanyarray(stars_gaia['DE_ICRS'])
dist_hyg = np.asanyarray(stars_gaia['Dist50']*1e3) # in parsecs
x_GCNS = np.asanyarray(stars_gaia['xcoord50']) # in parsecs 
y_GCNS = np.asanyarray(stars_gaia['ycoord50']) # in parsecs 
z_GCNS = np.asanyarray(stars_gaia['zcoord50']) # in parsecs 
print(f"Number of rows in the CSV: {len(stars_gaia)}")

# Filter for stars visible by the Green Bank Telescope
visible_stars = stars_gaia[(stars_gaia['Dist50']<=100000) & (stars_gaia['DE_ICRS'] >= -46) 
                           & (stars_gaia['DE_ICRS'] <= 90) & (stars_gaia['RA_ICRS'] >= 0) & (stars_gaia['RA_ICRS'] <=24)]

# Converting the ra and dec to radians and the distance to light years
ra_radians = np.radians(visible_stars['RA_ICRS'])  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_stars['DE_ICRS'])
stars_dist = visible_stars['Dist50']*3.26156

# Convert star positions to light-years
x_GCNS_ly = visible_stars['xcoord50']*3.26156 # in lightyears 
y_GCNS_ly = visible_stars['ycoord50']*3.26156
z_GCNS_ly = visible_stars['zcoord50']*3.26156
 #%%
 
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
y_sn = dist_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
z_sn = dist_sn *np.sin(np.radians(dec_sn))


event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)


# Calculate the distance between each star and the supernova
d2x = x_sn - x_GCNS_ly
d2y = y_sn - y_GCNS_ly
d2z = z_sn - z_GCNS_ly

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)


earth = 0

d3= stars_dist-earth

d1 = event_dist-earth

signal_time = (d2+d3)-d1

schedule = 1987+signal_time

visible_stars['schedule'] = 1987+round(signal_time)
#%%
visible_stars['Dist50'] = visible_stars['Dist50']*3.26156  # in light years

target_year = 2022
range_years = 10

# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year

# Filtering the catalogue bed on the target year
visible_stars = visible_stars[(visible_stars['schedule'] >= target_year - range_years) & 
                              (visible_stars['schedule'] <= target_year + range_years)]


# Correct coordinates for targets based on the year of observation

# coverting proper motion to correct coordinate
proper_ra_star = visible_stars['RA_ICRS'] * years_since / 3600  # Convert arcseconds to degrees
proper_dec_star = visible_stars['DE_ICRS'] * years_since / 3600  

# Calculate ra and dec using proepr motion
new_ra_star =visible_stars['RA_ICRS'] + proper_ra_star
new_dec_star = visible_stars['DE_ICRS'] + proper_dec_star

new_ra_rad_star = np.radians(new_ra_star)  
new_dec_rad_star = np.radians(new_dec_star)

# new star position
new_star_x = visible_stars['Dist50'] * np.cos(new_dec_rad_star) * np.cos(new_ra_rad_star)
new_star_y = visible_stars['Dist50'] * np.cos(new_dec_rad_star) * np.sin(new_ra_rad_star)
new_star_z = visible_stars['Dist50'] * np.sin(new_dec_rad_star)


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
c_sch = visible_stars.loc[visible_stars['schedule'] == visible_stars['schedule'].max()]
c_sch1  = visible_stars.loc[visible_stars['schedule'] == target_year]
c_sch2 = visible_stars.loc[visible_stars['schedule'] == visible_stars['schedule'].min()]

# ra and dec in terms of proper motion
proper_ra = float(c_sch['RA_ICRS'].values[0]) * years_since / 3600  # Convert arcseconds to degrees
proper_dec = float(c_sch['DE_ICRS'].values[0]) * years_since / 3600

new_ra = float(c_sch['RA_ICRS'].values[0]) + proper_ra
new_dec = float(c_sch['DE_ICRS'].values[0]) + proper_dec

new_ra_rad = np.radians(new_ra * 15)  
new_dec_rad = np.radians(new_dec)

# new position of star for ellipse
new_star_xe = float(c_sch['Dist50'].values[0]) * np.cos(new_dec_rad) * np.cos(new_ra_rad)
new_star_ye = float(c_sch['Dist50'].values[0]) * np.cos(new_dec_rad) * np.sin(new_ra_rad)
new_star_ze = float(c_sch['Dist50'].values[0]) * np.sin(new_dec_rad)

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

#%%

# Converting the coordinates of the data to light years
x_GCNS_ly_original = x_GCNS * 3.26156  # Convert parsecs to light-years
y_GCNS_ly_original = y_GCNS * 3.26156
z_GCNS_ly_original = z_GCNS * 3.26156
# Obtaining all the stars that are inside the ellipse can also use this code to obtain the stars outside the ellipse as well
inside_ellipse = (x_GCNS_ly_original**2 / a**2 + y_GCNS_ly_original**2 / b**2) <= 1

# Plot ellipe ontop of the ellipsoid 
plt.plot(x, y, color = 'red')
#plt.scatter(new_star_x, new_star_y, color = 'black', s=1)
plt.plot(sn_foci_x, sn_foci_y, 'bo')
plt.plot(e_foci_x, e_foci_y, 'go')
plt.scatter(x_GCNS_ly_original[inside_ellipse], y_GCNS_ly_original[inside_ellipse], s=1, color='pink', label="Stars") # plotting stars inside the ellipse
plt.scatter(x_GCNS_ly_original[~inside_ellipse], y_GCNS_ly_original[~inside_ellipse], s=1, color='orange', label="Stars")  # plotting stars outside the ellipse
plt.xlabel('lightyears')
plt.ylabel('lightyears')
plt.legend(['2022', 'earth','event', 'targets within the ellipse','targets outside of ellipse' ])
plt.axis('equal')
plt.xlim(-600, 600)
plt.ylim(-500, 400)
plt.show()