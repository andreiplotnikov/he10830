from astropy import units as u
from sunpy.net import attrs as a
import astropy.io.fits as fits
import astropy
import sunpy
import numpy as np
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import os
import matplotlib.pyplot as plt
import skimage
import scipy.stats
import scipy.ndimage
import drms
from sunpy.net import Fido
from reproject import reproject_interp
import sunpy.map
import datetime
import time

def limb_poly(x, ac):
    a = ac/np.sum(ac)
    s = np.zeros_like(x)
    mu = np.cos(x * np.pi * 0.5)
    for i in range(len(a)):
        s = s + a[i]*np.power(mu, i)
    return s
            

files_list = []

directory = 'C:\\results\\helium-ready\\'
for subd in os.listdir(directory):
    for subf in os.listdir(directory + subd):
        if 'fits' in subf:
            files_list.append(directory + subd + '\\' + subf)
            
small_list = []
big_list = []
small_files = []
big_files = []

for file in files_list:
    fixed_map = sunpy.map.Map(file)
    
    hpc_coords = sunpy.map.maputils.all_coordinates_from_map(fixed_map)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / fixed_map.rsun_obs
    r = r.value
    s = fixed_map.data
    
    s0 = np.median(s[r < 0.01])
    
    
    if (fixed_map.data.shape[0] < 300):
        small_list.append(s0)
        small_files.append(file)
        
    if (fixed_map.data.shape[0] > 600):
        big_list.append(s0)        
        big_files.append(file)
               

