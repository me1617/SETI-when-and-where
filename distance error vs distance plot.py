# -*- coding: utf-8 -*-
"""
Created on Fri Apr 25 16:23:22 2025

@author: meg16
"""



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
import setisynch.gaia as ss
import setisynch.hygdata_v3 as hs



votable = parse(r"C:\Users\meg16\Downloads\Data (6)\Data\GCNS_allcol.vot")  # Update with your actual file path
table = votable.get_first_table().to_table()

# Convert to Pandas DataFrame
gaia = table.to_pandas()
#print(gaia.columns)

# Filter for stars visible by the Green Bank Telescope
visible_stars = gaia[(gaia['Dist50']*3.26156e3<=100000) & (gaia['DE_ICRS'] >= -46) & (gaia['DE_ICRS'] <= 90) & (gaia['RA_ICRS'] >= 0) & (gaia['RA_ICRS'] <=24)]

ra_radians = np.radians(gaia['RA_ICRS'] * 15)  # Convert hours to degrees, then to radians
dec_radians = np.radians(gaia['DE_ICRS'])


planet_dist = gaia['Dist50']*3.26156e3


visible_planets_x = planet_dist * np.cos(dec_radians) * np.cos(ra_radians)
visible_planets_y = planet_dist  * np.cos(dec_radians) * np.sin(ra_radians)
visible_planets_z = planet_dist *np.sin(dec_radians)


ra_error = (gaia['e_RA_ICRS'])/gaia['RA_ICRS']
dec_error = (gaia['e_DE_ICRS'])/gaia['DE_ICRS']
#%%
dist_max_error = planet_dist*(1+(gaia['Dist84']-gaia['Dist50'])/gaia['Dist50'])
dist_min_error = planet_dist*(1+(gaia['Dist50']-gaia['Dist16'])/gaia['Dist50'])
import matplotlib.cm as cm
color_var = np.concatenate([ dist_min_error]) # You can also use `y_all` if you prefer
#plt.hist2d(planet_dist, dist_max_error, bins = 100)
#plt.hist2d(planet_dist, dist_min_error, bins = 100)
plt.scatter(planet_dist, dist_max_error, c=color_var, cmap=cm.plasma, s=10)
plt.scatter(planet_dist, dist_min_error, c=color_var, cmap=cm.plasma, s=10)

H, yedges, zedges = np.histogram2d(planet_dist, dist_max_error, bins = 100)
plt.pcolormesh(yedges, zedges, H)

