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
            
r_list = np.empty( (0,))
s_list = np.empty( (0,))
            
for file in files_list:
    fixed_map = sunpy.map.Map(file)
    
    if (fixed_map.data.shape[0] > 0):
        hpc_coords = sunpy.map.maputils.all_coordinates_from_map(fixed_map)
        r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / fixed_map.rsun_obs
        r = r.value
        s = fixed_map.data
        
        #center intensity
        s0 = np.median(s[r < 0.01])
        
        s = s/s0
        r = r[s > 0.6]
        s = s[s > 0.6]
        
        r_list = np.concatenate( (r_list, r))
        s_list = np.concatenate( (s_list, s))
        
print(r_list.shape)

points = np.vstack((r_list, s_list))

samples = 1000000

random_points = points[:, np.random.randint(0, r_list.shape[0], samples)]

plt.figure()
plt.subplot(121)
plt.hist2d(r_list, s_list, bins = 1000)

plt.subplot(122)
plt.hist2d(random_points[0], random_points[1], bins = 1000)


p0 = [ 0.59043513,  1.36563951, -2.6384554 ,  2.78392684, -1.10154608]

p0 = [ 0.64434808,  1.29852297, -2.62043326,  2.72090958, -1.04334737]

start = time.time()

k = scipy.optimize.basinhopping(lambda x: np.sum(np.power(limb_poly(random_points[0], x) - random_points[1], 2)), x0 = p0)

print('Time: ', time.time() - start)

print(k['fun'])

plt.hist2d(r_list, s_list, bins = 1000)
plt.plot(np.linspace(0, 1), limb_poly(np.linspace(0, 1), k['x']))

print(k['x']/np.sum(k['x']))