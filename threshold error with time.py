# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 17:00:21 2025

@author: meg16
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse
import setisynch.gaia as ss


on_01 = pd.read_csv(r"C:\Users\meg16\Downloads\Data (6)\Data\on_el_th0.1lyr.csv")
on_02 = pd.read_csv(r"C:\Users\meg16\Downloads\Data (3)\Data\on_el_th0.1lyr_2.csv")
on_03 = pd.read_csv(r"C:\Users\meg16\Downloads\Data (3)\Data\on_el_th0.1lyr_3.csv")
on_04 =  pd.read_csv(r"C:\Users\meg16\Downloads\Data (3)\Data\on_el_th0.1lyr_4.csv")

def error_data(on, SN1987A, target_year = 2096):
    # Filter for stars visible by the Green Bank Telescope
    visible_stars = on[(on['Dist50']*3.26156e3<=100000) & (on['DE_ICRS'] >= -46) & (on['DE_ICRS'] <= 90) & (on['RA_ICRS'] >= 0) & (on['RA_ICRS'] <=24)]
    
    ra_radians = np.radians(visible_stars['RA_ICRS'] * 15)  # Convert hours to degrees, then to radians
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
    
    
    
    planet_ra_new = visible_stars['RA_ICRS']*15
    planet_dec_new = visible_stars['DE_ICRS']
    
    ra_error = (visible_stars['e_RA_ICRS']- visible_stars['RA_ICRS'])/visible_stars['RA_ICRS']
    dec_error = (visible_stars['e_DE_ICRS'] - visible_stars['DE_ICRS'])/visible_stars['DE_ICRS']
    
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
    min_dist = (planet_dist * (1 + dist_min_error))
    max_dist = (planet_dist * (1 + dist_max_error))
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
    
    event_dist = np.sqrt(x_sn**2+y_sn**2)
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
    
    schedule_min = 1987+(signal_time*(1+min_signal_time))
    schedule_max = 1987+(signal_time*(1+max_signal_time))
    return schedule_min, schedule_max, min_dist, max_dist

SN1987A = ss.Event(ra=ss.hms2deg(5, 35, 28.03),
                      dec=ss.deg_arcmin_arcsec2deg(-69, 17, 11.79),
                      dist=ss.pc2ly(51.4e3),
                    event_year=ss.date2years('1987-02-24'),
                      event_name='SN1987A')

on_01_error = error_data(on_01, SN1987A)
on_02_error = error_data(on_02, SN1987A)
on_03_error = error_data(on_03, SN1987A)
on_04_error = error_data(on_04, SN1987A)

on_01_error_dist_min = on_01_error[2]
# Want to plot threshold, number of stars, schedule min
plt.hist(on_03_error[0], bins = 200, color  = 'green')
plt.hist(on_04_error[0], bins = 200, color= 'orange')
plt.hist(on_02_error[0], bins = 200, color = 'red')
plt.hist(on_01_error[0], bins = 200, color = 'blue')
plt.xlabel('Minimum Error Time')
plt.ylabel('Number of Stars')
import matplotlib.ticker as ticker

plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:.0f}'))
#%%
# Want to plot threshold, number of stars, schedule min
plt.plot(on_03_error[0],on_03_error[2],  color  = 'green')
#plt.plot(on_04_error[0], on_04_error[2],   color= 'orange')
#plt.plot(on_02_error[0],on_02_error[2] ,  color = 'red')
#plt.plot(on_01_error[0],on_01_error[2] , color = 'blue')
plt.xlabel('Minimum Error Time')
plt.ylabel('Minimum Distance Error')




    
#%%

    # # Create subplots~
    fig, ax = plt.subplots(3, 2, figsize=(12, 5))
    
    # Min error plot
    ax[0,0].plot(min_signal_time, dist_min_error,'.',   label='Min Error')
    ax[0,0].set_xlabel('Min Error Time')
    ax[0,0].set_ylabel('Min Error Distance')
    ax[0,0].legend()
    
    
    # Max error plot
    ax[0,1].plot(max_signal_time, dist_max_error,'.', label='Max Error', color='red')
    ax[0,1].set_xlabel('Max Error Time')
    ax[0,1].set_ylabel('Max Error Distance')
    ax[0,1].legend()
    
    ax[1,0].plot(max_signal_time, ra_error,'.', label='Max Error', color='red')
    ax[1,0].set_xlabel('Max Error Time')
    ax[1,0].set_ylabel('Error RA')
    ax[1,0].legend()
    
    ax[1,1].plot(min_signal_time, ra_error,'.', label='Min Error', color='red')
    ax[1,1].set_xlabel('Min Error Time')
    ax[1,1].set_ylabel('Error RA')
    ax[1,1].legend()
    
    ax[2,0].plot(max_signal_time, dec_error,'.', label='Max Error', color='red')
    ax[2,0].set_xlabel('Max Error Time')
    ax[2,0].set_ylabel('Error Dec')
    ax[2,0].legend()
    
    ax[2,1].plot(min_signal_time, dec_error,'.', label='Min Error', color='red')
    ax[2,1].set_xlabel('Min Error Time')
    ax[2,1].set_ylabel('Error Dec')
    ax[2,1].legend()
    
    
    
    
    plt.tight_layout()
    plt.show()
    min_schedule = 1987+(signal_time*(1+min_signal_time))
    max_schedule = 1987+(signal_time*(1+max_signal_time))
    
    visible_stars['min_schedule'] = min_schedule
    visible_stars['max_schedule'] = max_schedule
    
    
    #%%
    from datetime import datetime, timedelta
    
    def decimal_year_to_date(decimal_year):
        year = int(decimal_year)
        remainder = decimal_year - year
        start_of_year = datetime(year, 1, 1)
        days_in_year = (datetime(year + 1, 1, 1) - start_of_year).days
        # min_date = start_of_year + timedelta(days=remainder * days_in_year)
        # max_date = start_of_year + timedelta(days=remainder * days_in_year)
        date = start_of_year + timedelta(days=remainder * days_in_year)
        return date
    
    #%%
    
