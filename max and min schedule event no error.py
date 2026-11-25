# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 12:41:33 2025

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


event_dist = np.sqrt(x_sn**2+y_sn**2+z_sn**2)

d2x = x_sn-visible_planets_x
d2y  = y_sn-visible_planets_y
d2z = z_sn-visible_planets_z

d2 = np.sqrt(d2x**2+d2y**2+d2z**2)

earth = 0

d3 = planet_dist-earth

d1 = event_dist-earth

signal_time = (d2+d3)-d1
schedule = 1987+(signal_time)

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
plt.legend([target_year, 'stars', 'Event', 'Earth' ], fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.xlabel('lightyears', fontsize = 20)
plt.ylabel('lightyears', fontsize = 20)
plt.axis('equal')
plt.show()

#%%

#Calculation for the errors in ra, dec and distance
star_ra_new = np.radians(visible_stars['RA_ICRS'])
star_dec_new = np.radians(visible_stars['DE_ICRS'])

star_dist = visible_stars['Dist50']

#Calculating the percentage errors for ra, dec and distance
ra_error = ((visible_stars['e_RA_ICRS'])/3.6e6)/(visible_stars['RA_ICRS'])  # convert ra error from milliarcseconds to degrees
dec_error = ((visible_stars['e_DE_ICRS'])/3.6e6)/visible_stars['DE_ICRS']
dist_max_error = (visible_stars['Dist84']-visible_stars['Dist50'])/visible_stars['Dist50']
dist_min_error = (visible_stars['Dist50']-visible_stars['Dist16'])/visible_stars['Dist50']

# Applying the percentage error to the actual values of ra and dec
new_ra_error = star_ra_new*(1+ra_error)
new_dec_error = star_dec_new*(1+dec_error)


x_sn = dist_sn * np.cos(np.radians(star_dec_new)) * np.cos(np.radians(star_ra_new))
y_sn = dist_sn  * np.cos(np.radians(star_dec_new)) * np.sin(np.radians(star_ra_new))
z_sn = dist_sn *np.sin(np.radians(star_dec_new))

visible_stars_x = star_dist * np.cos(np.radians(star_dec_new)) * np.cos(np.radians(star_ra_new))
visible_stars_y = star_dist * np.cos(np.radians(star_dec_new)) * np.sin(np.radians(star_ra_new))
visible_stars_z = star_dist * np.sin(np.radians(star_dec_new))


def error_coordinates(dist_error ,ra, dec, ra_e, dec_e):
    cos_r_error = np.sin(np.radians(ra))*np.radians(ra_e) # cosine error propagating equation
    cos_r_actual = np.cos(np.radians(ra))
    cos_d_error = np.sin(np.radians(dec))*np.radians(dec_e) # cosine error propagating equation
    cos_d_actual = np.cos(np.radians(dec))
    sin_r_error = np.cos(np.radians(ra))*np.radians(ra_e) # sine error propagating equation
    sin_r_actual = np.sin(np.radians(ra))
    sin_d_error = np.cos(np.radians(dec))*np.radians(dec_e) # sine error propagating equation
    sin_d_actual = np.sin(np.radians(dec))
    x_error = np.sqrt((dist_error)**2 + (cos_d_error/cos_d_actual)**2 + (cos_r_error/cos_r_actual)**2)
    y_error = np.sqrt((dist_error)**2 + (cos_d_error/cos_d_actual)**2 + (sin_r_error/sin_r_actual)**2)
    z_error = np.sqrt((dist_error)**2 + (sin_r_error/sin_r_actual)**2)
    return x_error, y_error, z_error
    

def distance(x_1, y_1, x_2, y_2):
    return np.sqrt((x_1-x_2)**2+(y_1-y_2)**2)

dist_max = (star_dist * (1 + dist_max_error))
dist_min = (star_dist * (1 + dist_min_error))

max_stars_x,max_stars_y,max_stars_z = error_coordinates(dist_max_error,star_ra_new, star_dec_new, ra_error, dec_error)
min_stars_x, min_stars_y ,min_stars_z  = error_coordinates(dist_min_error,star_ra_new, star_dec_new, ra_error, dec_error)

d2 = distance(x_sn, y_sn, visible_stars_x, visible_stars_y)

# Calculating the error maximum distance between the event and the stars
max_d2 = (distance(x_sn, y_sn, max_stars_x*visible_stars_x, max_stars_y*visible_stars_y)-d2)/d2

# Calculating the minimum percentage error distance between the event and the stars
min_d2 = (distance(x_sn, y_sn, min_stars_x*visible_stars_x, min_stars_y*visible_stars_y)-d2)/d2

event_dist = np.sqrt(x_sn**2+y_sn**2)
max_event_dist = np.sqrt(x_sn**2+y_sn**2)
min_event_dist = np.sqrt(x_sn**2+y_sn**2)

earth = 0
d3 = star_dist-earth
d1 = event_dist-earth

# Calculating the percentage error distance between the stars and earth
max_d3 = (((star_dist * (1 + dist_max_error))-earth)-d3)/d3
min_d3 = (d3-((star_dist * (1 + dist_min_error))-earth))/d3

# Calculating the percentage error distance between the event and earth
max_d1 = ((max_event_dist) - d1) / d1
min_d1 =  (d1-min_event_dist) / d1

signal_time = (d2+d3)-d1
# Calculating the signal_time error using prop
min_signal_time =np.sqrt(min_d2**2+min_d3**2+min_d1**2)#/signal_time
max_signal_time =np.sqrt(max_d2**2+max_d3**2+max_d1**2)#/signal_time


# Create subplots
fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# Min and max error plot time and distance
ax[0].plot(min_signal_time*100, dist_min_error*100,'.')
ax[0].set_xlabel('Min Error Time', fontsize = 20)
ax[0].set_ylabel('Min Error Distance', fontsize = 20)
ax[0].tick_params(axis='x', labelsize=20)
ax[0].tick_params(axis='y', labelsize=20)

ax[1].plot(max_signal_time*100, dist_max_error*100,'.', label='Max Error', color='red')
ax[1].set_xlabel('Max Error Time', fontsize = 20)
ax[1].set_ylabel('Max Error Distance', fontsize = 20)
ax[1].tick_params(axis='x', labelsize=20)
ax[1].tick_params(axis='y', labelsize=20)

plt.tight_layout()
plt.show()

# Max and min error plot time and ra
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(max_signal_time*100, ra_error*100,'.', label='Max Error', color='red')
ax[0].set_xlabel('Max Error Time', fontsize = 20)
ax[0].set_ylabel('Error RA', fontsize = 20)
ax[0].tick_params(axis='x', labelsize=20)
ax[0].tick_params(axis='y', labelsize=20)

ax[1].plot(min_signal_time*100, ra_error*100,'.', label='Min Error', color='red')
ax[1].set_xlabel('Min Error Time', fontsize = 20)
ax[1].set_ylabel('Error RA', fontsize = 20)
ax[1].tick_params(axis='x', labelsize=20)
ax[1].tick_params(axis='y', labelsize=20)

plt.tight_layout()
plt.show()

# Max and min error plot time and dec
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(max_signal_time*100, dec_error*100,'.', label='Max Error', color='red')
ax[0].set_xlabel('Max Error Time', fontsize = 20)
ax[0].set_ylabel('Error Dec', fontsize = 20)
ax[0].tick_params(axis='x', labelsize=20)
ax[0].tick_params(axis='y', labelsize=20)

ax[1].plot(min_signal_time*100, dec_error*100,'.', label='Min Error', color='red')
ax[1].set_xlabel('Min Error Time', fontsize = 20)
ax[1].set_ylabel('Error Dec', fontsize = 20)
ax[1].tick_params(axis='x', labelsize=20)
ax[1].tick_params(axis='y', labelsize=20)


plt.tight_layout()
plt.show()

# Creating the minimum and maximum schedule
min_signal_error = (signal_time*(1+min_signal_time))
min_schedule = schedule-min_signal_error
max_signal_error = (signal_time*(1+max_signal_time))
max_schedule = schedule+max_signal_error

# Adding the schedules to the dataset
visible_stars['min_schedule'] = min_schedule
visible_stars['max_schedule'] = max_schedule


#%%
# Filtering the catalogue bed on the target year

target_year_min = 2025
range_years = 10
years_since_min = target_year_min - SN_year
# Filtering the catalogue on the target year for the maximum schedule
visible_stars_min = visible_stars[(visible_stars['min_schedule'] >= target_year_min - range_years) & 
                              (visible_stars['min_schedule'] <= target_year_min + range_years)]

# Obtain and recalculate the ra, dec and distance values and their errors
star_ra_new_min = np.radians(visible_stars_min['RA_ICRS'])
star_dec_new_min = np.radians(visible_stars_min['DE_ICRS'])
ra_error_min = ((visible_stars_min['e_RA_ICRS'])/3.6e6)/(visible_stars_min['RA_ICRS'])  # convert ra error from milliarcseconds to degrees
dec_error_min = ((visible_stars_min['e_DE_ICRS'])/3.6e6)/visible_stars_min['DE_ICRS']
dist_min_error = (visible_stars_min['Dist50']-visible_stars_min['Dist16'])/visible_stars_min['Dist50']
dist_max = visible_stars_min['Dist50']*(1+dist_min_error)

# Applying the time to the ra and dec
ra_star_year_min = (visible_stars_min['RA_ICRS']*(1+ra_error_min)) * years_since_min / 3600  # Convert arcseconds to degrees
dec_star_year_min = (visible_stars_min['DE_ICRS']*(1+dec_error_min)) * years_since_min / 3600  
new_ra_star_min = np.radians(ra_star_year_min)  
new_dec_star_min = np.radians(dec_star_year_min)

# new star position
new_star_x_min, new_star_y_min, new_star_z_min = error_coordinates(dist_min,star_ra_new_min, star_dec_new_min, new_ra_star_min,new_dec_star_min)

# starting to create the ellipse
target_schedules_min = target_year_min

# Setting up for ellipse
# Correct coordinates for targets based on the year of observation for plotting the ellipse
c_sch1_min = visible_stars_min.loc[(visible_stars_min['min_schedule'].astype(int)) == (visible_stars_min['min_schedule'].astype(int)).max()]
c_sch_min  = visible_stars_min.loc[(visible_stars_min['min_schedule'].astype(int)) == target_year_min]
c_sch2_min = visible_stars_min.loc[(visible_stars_min['min_schedule'].astype(int)) == (visible_stars_min['min_schedule'].astype(int)).min()]

# Calculating the ra, dec and distance error for c_sch
c_sch_ra_min = (c_sch_min['e_RA_ICRS']/3.6e6)/c_sch_min['RA_ICRS']
c_sch_dec_min = (c_sch_min['e_DE_ICRS']/3.6e6)/c_sch_min['DE_ICRS']
c_sch_min_dist = c_sch_min['Dist50']*(1+((c_sch_min['Dist50']-c_sch_min['Dist16'])/c_sch_min['Dist50']))
ra_year_min = (c_sch_min['RA_ICRS']*(1+c_sch_ra_min)) * years_since_min / 3600  # Convert arcseconds to degrees
dec_year_min = (c_sch_min['DE_ICRS']*(1+c_sch_dec_min))* years_since_min / 3600
new_ra_min = np.radians(ra_year_min)  
new_dec_min = np.radians(dec_year_min)

# new position of star for ellipse
new_star_xe_min,new_star_ye_min,new_star_ze_min =  error_coordinates(c_sch_min_dist,np.radians(c_sch_min['RA_ICRS']),np.radians(c_sch_min['DE_ICRS']) ,
                                                                     new_ra_min,new_dec_min)

# Calculate distances from the event
event_2_starx_min = x_sn - float(new_star_xe_min.values[0])
event_2_stary_min = y_sn - float(new_star_ye_min.values[0])
event_2_starz_min = z_sn - float(new_star_ze_min.values[0])
star_2_earth_min = float(c_sch_min_dist.values[0])

event_2_star_min = np.sqrt((float(new_star_xe_min.values[0]) - x_sn)**2 + (float(new_star_ye_min.values[0]) - y_sn)**2)

# foci
sn_foci_x_min,sn_foci_y_min  = x_sn, y_sn
e_foci_x, e_foci_y = 0, 0


x0_min = (sn_foci_x_min + e_foci_x) / 2                       # Center x-value
y0_min = (sn_foci_y_min + e_foci_y) / 2                       # Center y-value
f_min = np.sqrt((sn_foci_x - x0)**2 + (sn_foci_y - y0)**2) # Distance from center to focus

c_min = event_2_star+star_2_earth_min

# Compute ellipse parameters
a_min = c_min / 2                                # Semimajor axis
b_min = np.sqrt(a_min**2 - f_min**2)                 # Semiminor axis
phi_min = np.arctan2((e_foci_y - sn_foci_y_min), (e_foci_x - sn_foci_x_min))   # Angle betw major y-axis and x-axis

# Parametric plot in t
resolution = 10000
t = np.linspace(0, 2*np.pi, resolution)
x_min = x0 + a_min * np.cos(t) * np.cos(phi) - b_min * np.sin(t) * np.sin(phi)
y_min = y0 + a_min * np.cos(t) * np.sin(phi) + b_min * np.sin(t) * np.cos(phi)


fig = plt.figure()
# Plot ellipse
plt.plot(x_min, y_min, color = 'red')

plt.scatter(new_star_xe_min, new_star_ye_min, color = 'purple')
# Show focii
plt.plot(sn_foci_x, sn_foci_y, 'bo')
plt.plot(e_foci_x, e_foci_y, 'go')
plt.legend([target_year_min, 'stars', 'Event', 'Earth' ], fontsize = '20')
plt.xlabel('lightyears')
plt.ylabel('lightyears')
plt.axis('equal')
plt.show()
#%%

target_year_max = 2025
range_years = 10
years_since_max = target_year_max - SN_year
# Filtering the catalogue on the target year for the maximum schedule
visible_stars_max = visible_stars[(visible_stars['max_schedule'] >= target_year_max - range_years) & 
                              (visible_stars['max_schedule'] <= target_year_max + range_years)]

# Obtain and recalculate the ra, dec and distance values and their values
star_ra_new_max = np.radians(visible_stars_max['RA_ICRS'])
star_dec_new_max = np.radians(visible_stars_max['DE_ICRS'])
ra_error_max = (visible_stars_max['e_RA_ICRS']/3.6e6)/visible_stars_max['RA_ICRS']
dec_error_max = (visible_stars_max['e_DE_ICRS']/3.6e6)/visible_stars_max['DE_ICRS']
dist_max_error = (visible_stars_max['Dist84']-visible_stars_max['Dist50'])/visible_stars_max['Dist50']
dist_max = visible_stars_max['Dist50']*(1+dist_max_error)

# Applying the time to the ra and dec
ra_star_year_max = (visible_stars_max['RA_ICRS']*(1+ra_error_max))  * years_since_max / 3600  # Convert arcseconds to degrees
dec_star_year_max = (visible_stars_max['DE_ICRS']*(1+dec_error_max))  * years_since_max / 3600  
new_ra_star_max = np.radians(ra_star_year_max)  
new_dec_star_max = np.radians(dec_star_year_max)

# Calculating the new stars max error coordinates
new_star_x_max, new_star_y_max, new_star_z_max = error_coordinates(dist_max,star_ra_new_max, star_dec_new_max, new_ra_star_max,new_dec_star_max)

# Setting up for ellipse

# Correct coordinates for targets based on the year of observation for plotting the ellipse
c_sch1_max = visible_stars_max.loc[(visible_stars_max['max_schedule'].astype(int)) == (visible_stars_max['max_schedule'].astype(int)).max()]
c_sch_max  = visible_stars_max.loc[(visible_stars_max['max_schedule'].astype(int)) == target_year_max]
c_sch2_max = visible_stars_max.loc[(visible_stars['max_schedule'].astype(int)) == (visible_stars_max['max_schedule'].astype(int)).min()]

# Calculating the ra, dec and distance errors for the new filtered stars for the year of interest
c_sch_ra_max = (c_sch_max['e_RA_ICRS']/3.6e6)/c_sch_max['RA_ICRS']
c_sch_dec_max = (c_sch_max['e_DE_ICRS']/3.6e6)/c_sch_max['DE_ICRS']

c_sch_max_dist = c_sch_max['Dist50']*(1+(c_sch_max['Dist84']-c_sch_max['Dist50'])/c_sch_max['Dist50'])

ra_year_max = (c_sch_max['RA_ICRS']*(1+c_sch_ra_max))  * years_since_max / 3600  # Convert arcseconds to degrees
dec_year_max = (c_sch_max['RA_ICRS']*(1+c_sch_dec_max))* years_since_max / 3600

new_ra_max = np.radians(ra_year_max)  
new_dec_max = np.radians(dec_year_max)

# New position of star for ellipse
new_star_xe_max = c_sch_max_dist* np.cos(new_dec_max) * np.cos(new_ra_max)
new_star_ye_max = c_sch_max_dist * np.cos(new_dec_max) * np.sin(new_ra_max)
new_star_ze_max = c_sch_max_dist * np.sin(new_dec_max)

# Calculate distances from the event to stars and stars to earth
event_2_starx_max = x_sn - float(new_star_xe_max.values[0])
event_2_stary_max = y_sn - float(new_star_ye_max.values[0])
event_2_starz_max = z_sn - float(new_star_ze_max.values[0])
star_2_earth_max = float(c_sch_max_dist.values[0])

event_2_star_max = np.sqrt((float(new_star_xe_max.values[0]) - x_sn)**2 + (float(new_star_ye_max.values[0]) - y_sn)**2)

# Foci
sn_foci_x_max,sn_foci_y_max  = x_sn, y_sn
e_foci_x, e_foci_y = 0, 0


# Calculate ellipse parameters
x0_max = (sn_foci_x_max + e_foci_x) / 2                       # Center x-value
y0_max = (sn_foci_y_max + e_foci_y) / 2                       # Center y-value
f_max = np.sqrt((sn_foci_x_max - x0_max)**2 + (sn_foci_y_max - y0_max)**2) # Distance from center to focus

c_max = event_2_star+star_2_earth_max

a_max = c_max / 2                                # Semimajor axis
b_max= np.sqrt(a_max**2 - f**2)                 # Semiminor axis
phi_max = np.arctan2((e_foci_y - sn_foci_y_max), (e_foci_x - sn_foci_x_max))   # Angle betw major y-axis and x-axis

# Parametric plot in t
resolution = 10000
t = np.linspace(0, 2*np.pi, resolution)

#Calculating the coordinates for the ellipse
x_max = x0 + a_max * np.cos(t) * np.cos(phi) - b_max * np.sin(t) * np.sin(phi)
y_max = y0 + a_max * np.cos(t) * np.sin(phi) + b_max * np.sin(t) * np.cos(phi)


fig = plt.figure()
# Plot ellipse
plt.plot(x,y, color = 'red' )
plt.plot(x_max, y_max, color = 'orange')
plt.plot(x_min, y_min, color = 'fuchsia')
plt.scatter(new_star_xe, new_star_ye, color = 'purple')
plt.scatter(new_star_xe_max, new_star_ye_max, color = 'magenta')
plt.scatter(new_star_xe_min, new_star_ye_min, color = 'violet')

# Show focii
plt.plot(sn_foci_x, sn_foci_y, color = 'blue', marker  = 'x')
plt.plot(sn_foci_x, sn_foci_y, color = 'cyan', marker  = 'x' )
plt.plot(sn_foci_x, sn_foci_y, color = 'indigo', marker  = 'x')


plt.plot(e_foci_x, e_foci_y, 'go')
plt.legend([2025, 'Max Error 2025','Min Error 2025',  'stars','Max Error Stars' ,'Min Error Stars', 'Event', 'Max Error Event', 'Min Error Event',  'Earth' ], fontsize = '20')
plt.xticks(fontsize = '20')
plt.yticks(fontsize = '20')
plt.xlabel('X (ly)', fontsize = '20')
plt.ylabel('Y (ly)', fontsize = '20')
plt.axis('equal')
plt.show()



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
min_date = (c_sch_min['min_schedule']).apply(real_date)
max_date = (c_sch_max['max_schedule']).apply(real_date)
date = (c_sch_max['schedule']).apply(real_date)

# Then saving them in the file 
c_sch_max['date'] = date
c_sch_min['min_date'] = min_date
c_sch_max['max_date'] = max_date
# Defining the dile data in a zip to then create the table, with only the ra, dec, distance and date values including their errors
final_data = list(zip( round(c_sch_max[:10]['RA_ICRS'], 3),round(c_sch_max[:10]['e_RA_ICRS']/3.6, 3), round(c_sch_max[:10]['DE_ICRS'], 3), round(c_sch_max[:10]['e_DE_ICRS'], 3),round(c_sch_max[:10]['Dist16'], 3),
                      round(c_sch_max[:10]['Dist50'], 3),round(c_sch_max[:10]['Dist84'], 3),round(c_sch_max[:10]['date'], 3), round(c_sch_min[:10]['min_date'], 3), round(c_sch_max[:10]['max_date'], 3)))
# Table headers
headers = [ "RA (Degrees)","RA Error (Degrees)", "Dec (Degrees)", "Dec Error (Degrees)", "Min Distance (Lyrs)", "Distance (Lyrs)", "Max distance (Lyrs)", "Date", "Min Date", "Max Date"]


# Create a figure and axis
import matplotlib.pyplot as plt

# Create figure and axis
fig, ax = plt.subplots(figsize=(30, 35))  # Large figure
ax.set_axis_off()

# Create the table
table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')

# Adjust font size for all cells
table.auto_set_font_size(False)
table.set_fontsize(15)  # Set your desired font size

# Optionally adjust row/column sizes
table.scale(1.2, 2)  # (col_width, row_height)
plt.tight_layout()
plt.show()

#%%

# Extracting the first 30 rows for the c_sch_min calculation and applying the time function to the original and minimum schedule
min_date = (c_sch_min['min_schedule']).apply(real_date)
date = (c_sch_min['schedule']).apply(real_date)

# Then saving them in the file 
c_sch_min['date'] = date
c_sch_min['min_date'] = min_date

# Defining the dile data in a zip to then create the table, with only the ra, dec, distance and date values including their errors
final_data = list(zip( c_sch_min[:10]['RA_ICRS'],c_sch_min[:10]['e_RA_ICRS']/3.6, c_sch_min[:10]['DE_ICRS'], c_sch_min[:10]['e_DE_ICRS'],c_sch_min[:10]['Dist16'],
                      c_sch_min[:10]['Dist50'],c_sch_min[:10]['date'], c_sch_min[:10]['min_date'] ))
# Table headers
headers = [ "RA (Degrees)","RA Error (Degrees)", "Dec (Degrees)", "Dec Error (Degrees)", "Min Distance (ly)", "Distance (ly)", "Date", "Min Date"]


# Create a figure and axis
import matplotlib.pyplot as plt

# Create figure and axis
fig, ax = plt.subplots(figsize=(30, 35))  # Large figure
ax.set_axis_off()

# Create the table
table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')

# Adjust font size for all cells
table.auto_set_font_size(False)
table.set_fontsize(15)  # Set your desired font size

# Optionally adjust row/column sizes
table.scale(1.2, 2)  # (col_width, row_height)
plt.tight_layout()
plt.show()

#%%

# Extracting the first 30 rows for the c_sch_min calculation and applying the time function to the original and minimum schedule
max_date = (c_sch_max['max_schedule']).apply(real_date)
date = (c_sch_max['schedule']).apply(real_date)

# Then saving them in the file 
c_sch_max['date'] = date
c_sch_max['max_date'] = max_date
# Defining the dile data in a zip to then create the table, with only the ra, dec, distance and date values including their errors
final_data = list(zip( c_sch_max[:10]['RA_ICRS'],c_sch_max[:10]['e_RA_ICRS']/3.6, c_sch_max[:10]['DE_ICRS'], c_sch_max[:10]['e_DE_ICRS'],
                      c_sch_max[:10]['Dist50'],c_sch_max[:10]['Dist84'],c_sch_max[:10]['date'], c_sch_max[:10]['max_date'] ))
# Table headers
headers = [ "RA (Degrees)","RA Error (Degrees)", "Dec (Degrees)", "Dec Error (Degrees)", "Distance (Lyrs)", "Max distance (Lys)", "Date", "Max Date"]


# Create a figure and axis
import matplotlib.pyplot as plt

# Create figure and axis
fig, ax = plt.subplots(figsize=(30, 35))  # Large figure
ax.set_axis_off()

# Create the table
table = ax.table(cellText=final_data, colLabels=headers, cellLoc='center', loc='center')

# Adjust font size for all cells
table.auto_set_font_size(False)
table.set_fontsize(15)  # Set your desired font size

# Optionally adjust row/column sizes
table.scale(1.2, 2)  # (col_width, row_height)
plt.tight_layout()
plt.show()

# Calculating the maximum and minimum errors on the coordinates for the stars
def coordinates(dist, ra, dec):
    x_stars = dist * np.cos(np.radians(dec)) * np.cos(np.radians(ra))
    y_stars = dist * np.cos(np.radians(dec)) * np.sin(np.radians(ra))
    z_stars = dist * np.sin(np.radians(dec))
    return x_stars,y_stars,z_stars

