# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 16:08:21 2025

@author: meg16
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse
import setisynch.gaia as ss


# Parse the VOTable
votable = parse(r"C:\Users\meg16\Downloads\Data (6)\Data\GCNS_allcol.vot")  # Update with your actual file path
table = votable.get_first_table().to_table()
# Convert to Pandas DataFrame
gaia = table.to_pandas()
print(gaia.columns)


# Filter for stars visible by the Green Bank Telescope
visible_stars = gaia[(gaia['Dist50']*3.26156e3<=100000) & (gaia['DE_ICRS'] >= -46) & (gaia['DE_ICRS'] <= 90) & (gaia['RA_ICRS'] >= -360) & (gaia['RA_ICRS'] <=360)]

ra_radians = np.radians(visible_stars['RA_ICRS'])  # Convert hours to degrees, then to radians
dec_radians = np.radians(visible_stars['DE_ICRS'])


planet_dist = visible_stars['Dist50']*3.26156e3


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


event_dist = np.sqrt(x_sn**2+y_sn**2)

d2x = x_sn-visible_planets_x
d2y  = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z

d2 = np.sqrt(d2x**2+d2y**2)

earth = 0

d3 = planet_dist-earth

d1 = event_dist-earth

signal_time = (d2+d3)-d1
schedule = 1987+np.round(signal_time)

visible_stars['schedule'] = 1987+(signal_time)

visible_stars['signal_time'] = signal_time
visible_stars['Dist50'] = visible_stars['Dist50']*3.26156e3 # changing dataset so in light years instead of parsecs
visible_stars['Dist84'] = visible_stars['Dist84']*3.26156e3 
visible_stars['Dist16'] = visible_stars['Dist16']*3.26156e3 


target_year = 2025
range_years = 10

# year of event and when the targets will see event
SN_year = 1987  
years_since = target_year - SN_year

#%%
# Filtering the catalogue bed on the target year
visible_stars = visible_stars[(visible_stars['schedule'] >= target_year - range_years) & 
                              (visible_stars['schedule'] <= target_year + range_years)]


# Correct coordinates for targets based on the year of observation

# coverting proper motion to correct coordinate
ra_star_year = visible_stars['RA_ICRS'] * years_since / 3600  # Convert arcseconds to degrees
dec_star_year = visible_stars['DE_ICRS'] * years_since / 3600  

# Calculate ra and dec using proepr motion
new_ra_star =visible_stars['RA_ICRS'] + ra_star_year
new_dec_star = visible_stars['DE_ICRS'] + dec_star_year

new_ra_star = np.radians(new_ra_star)  
new_dec_star = np.radians(new_dec_star)

# new star position
new_star_x = visible_stars['Dist50'] * np.cos(new_dec_star) * np.cos(new_ra_star)
new_star_y = visible_stars['Dist50'] * np.cos(new_dec_star) * np.sin(new_ra_star)
new_star_z = visible_stars['Dist50'] * np.sin(new_dec_star)


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

# Correct coordinates for targets based on the year of observation for plotting the ellipse
c_sch = visible_stars.loc[(visible_stars['schedule'].astype(int)) == (visible_stars['schedule'].astype(int)).max()]
c_sch1  = visible_stars.loc[(visible_stars['schedule'].astype(int)) == target_year]
c_sch2 = visible_stars.loc[(visible_stars['schedule'].astype(int)) == (visible_stars['schedule'].astype(int)).min()]

# ra and dec in terms of proper motion
ra_year = c_sch['RA_ICRS'] * years_since / 3600  # Convert arcseconds to degrees
dec_year = c_sch['DE_ICRS']* years_since / 3600

new_ra = c_sch['RA_ICRS'] + ra_year
new_dec = c_sch['DE_ICRS'] + dec_year

new_ra = np.radians(new_ra)  
new_dec = np.radians(new_dec)

# new position of star for ellipse
new_star_xe = c_sch['Dist50'] * np.cos(new_dec) * np.cos(new_ra)
new_star_ye = c_sch['Dist50'] * np.cos(new_dec) * np.sin(new_ra)
new_star_ze = c_sch['Dist50'] * np.sin(new_dec)

# Calculate distances from the event
event_2_starx = x_sn - float(new_star_xe.values[0])
event_2_stary = y_sn - float(new_star_ye.values[0])
event_2_starz = z_sn - float(new_star_ze.values[0])
star_2_earth = float(c_sch['Dist50'].values[0])

event_2_star = np.sqrt((float(new_star_xe.values[0]) - x_sn)**2 + (float(new_star_ye.values[0]) - y_sn)**2)

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
#%%



planet_ra_new = np.radians(visible_stars['RA_ICRS'])
planet_dec_new = np.radians(visible_stars['DE_ICRS'])

planet_dist = visible_stars['Dist50']

ra_error = ((visible_stars['e_RA_ICRS'])/3.6)/(visible_stars['RA_ICRS'])
dec_error = ((visible_stars['e_DE_ICRS'])/3.6)/visible_stars['DE_ICRS']

dist_max_error = (visible_stars['Dist84']-visible_stars['Dist50'])/visible_stars['Dist50']
dist_min_error = (visible_stars['Dist50']-visible_stars['Dist16'])/visible_stars['Dist50']

new_ra_error = planet_ra_new*(1+ra_error)

new_dec_error = planet_ra_new*(1+dec_error)



x_sn = dist_sn * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
y_sn = dist_sn  * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
z_sn = dist_sn *np.sin(np.radians(planet_dec_new))

visible_planets_x = planet_dist * np.cos(np.radians(planet_dec_new)) * np.cos(np.radians(planet_ra_new))
visible_planets_y = planet_dist * np.cos(np.radians(planet_dec_new)) * np.sin(np.radians(planet_ra_new))
visible_planets_z = planet_dist * np.sin(np.radians(planet_dec_new))


max_planets_x = (planet_dist * (1 + dist_max_error)) * np.cos(np.radians(new_dec_error)) * np.cos(np.radians(new_ra_error))
max_planets_y = (planet_dist * (1 + dist_max_error)) * np.cos(np.radians(new_dec_error)) * np.sin(np.radians(new_ra_error))
max_planets_z = (planet_dist * (1 + dist_max_error)) * np.sin(np.radians(new_dec_error))

min_planets_x = (planet_dist * (1 + dist_min_error)) * np.cos(np.radians(new_dec_error)) * np.cos(np.radians(new_ra_error))
min_planets_y = (planet_dist * (1 + dist_min_error)) * np.cos(np.radians(new_dec_error)) * np.sin(np.radians(new_ra_error))
min_planets_z = (planet_dist * (1 + dist_min_error)) * np.sin(np.radians(new_dec_error))

max_sn = (dist_sn+(1.2e3*3.26156))
max_x_sn = max_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
max_y_sn = max_sn  * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
max_z_sn = max_sn * np.sin(np.radians(dec_sn))

max_d2x = max_x_sn - max_planets_x
max_d2y = max_y_sn - max_planets_y
max_d2z = max_z_sn - max_planets_z

d2x = x_sn-visible_planets_x
d2y = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z


d2 = np.sqrt(d2x**2+d2y**2)

max_d2 = (np.sqrt(max_d2x**2+max_d2y**2)-d2)/d2

min_sn = (dist_sn-(1.2e3*3.26156))

min_x_sn = min_sn * np.cos(np.radians(dec_sn)) * np.cos(np.radians(ra_sn))
min_y_sn = min_sn * np.cos(np.radians(dec_sn)) * np.sin(np.radians(ra_sn))
min_z_sn = min_sn * np.sin(np.radians(dec_sn))

min_d2x = min_x_sn - min_planets_x
min_d2y = min_y_sn - min_planets_y
min_d2z = min_z_sn - min_planets_z

min_d2 = (d2-np.sqrt(min_d2x**2+min_d2y**2))/d2

event_dist = np.sqrt(x_sn**2+y_sn**2)
max_event_dist = np.sqrt(max_x_sn**2+max_y_sn**2)
min_event_dist = np.sqrt(min_x_sn**2+min_y_sn**2)

earth = 0

d3 = planet_dist-earth

d1 = event_dist-earth


max_d3 = (((planet_dist * (1 + dist_max_error))+earth)-d3)/d3
min_d3 = (((planet_dist * (1 + dist_min_error))+earth)-d3)/d3


max_d1 = ((max_event_dist) - d1) / d1
min_d1 = (d1 - (min_event_dist)) / d1

signal_time = (d2+d3)-d1
min_signal_time = min_d2 + min_d3 - min_d1
max_signal_time = max_d2 + max_d3 - max_d1


# # Create subplots~
fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# Min error plot
ax[0].plot(min_signal_time, dist_min_error,'.')
ax[0].set_xlabel('Min Error Time', fontsize = 20)
ax[0].set_ylabel('Min Error Distance', fontsize = 20)
ax[0].tick_params(axis='x', labelsize=20)
ax[0].tick_params(axis='y', labelsize=20)


# Max error plot
ax[1].plot(max_signal_time, dist_max_error,'.', label='Max Error', color='red')
ax[1].set_xlabel('Max Error Time', fontsize = 20)
ax[1].set_ylabel('Max Error Distance', fontsize = 20)
ax[1].tick_params(axis='x', labelsize=20)
ax[1].tick_params(axis='y', labelsize=20)

plt.tight_layout()
plt.show()
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(max_signal_time, ra_error,'.', label='Max Error', color='red')
ax[0].set_xlabel('Max Error Time', fontsize = 20)
ax[0].set_ylabel('Error RA', fontsize = 20)
ax[0].tick_params(axis='x', labelsize=20)
ax[0].tick_params(axis='y', labelsize=20)

ax[1].plot(min_signal_time, ra_error,'.', label='Min Error', color='red')
ax[1].set_xlabel('Min Error Time', fontsize = 20)
ax[1].set_ylabel('Error RA', fontsize = 20)
ax[1].tick_params(axis='x', labelsize=20)
ax[1].tick_params(axis='y', labelsize=20)

plt.tight_layout()
plt.show()

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(max_signal_time, dec_error,'.', label='Max Error', color='red')
ax[0].set_xlabel('Max Error Time', fontsize = 20)
ax[0].set_ylabel('Error Dec', fontsize = 20)
ax[0].tick_params(axis='x', labelsize=20)
ax[0].tick_params(axis='y', labelsize=20)

ax[1].plot(min_signal_time, dec_error,'.', label='Min Error', color='red')
ax[1].set_xlabel('Min Error Time', fontsize = 20)
ax[1].set_ylabel('Error Dec', fontsize = 20)
ax[1].tick_params(axis='x', labelsize=20)
ax[1].tick_params(axis='y', labelsize=20)


plt.tight_layout()
plt.show()
min_signal_error = (signal_time*(1+min_signal_time))
min_schedule = 1987-min_signal_error
max_schedule = (signal_time*(1+max_signal_time))-1987

visible_stars['min_schedule'] = min_schedule
visible_stars['max_schedule'] = max_schedule