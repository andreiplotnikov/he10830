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
import shutil

directory = 'C:\\results\\helium\\'
check_directory = 'C:\\results\\helium-good\\'

save_path = 'C:\\results\\helium-ready\\'

files_list = []
files_names = []
files_copydir = []

for root, dirs, files in os.walk(directory, topdown = False):
    for name in files:
        files_names.append(name)
        files_copydir.append(os.path.basename(os.path.dirname(os.path.join(root, name))))
        files_list.append(os.path.join(root, name))
        
files_checklist = []

for root, dirs, files in os.walk(check_directory, topdown = False):
    for d in dirs:
        try:
            os.mkdir(save_path + '%s' % (d))
        except:
            pass
    for name in files:
        files_checklist.append(name)

       
for file_num in range(len(files_list)):
    if files_names[file_num] in files_checklist:
        shutil.copyfile(files_list[file_num], save_path + files_copydir[file_num] + '\\' + files_names[file_num])
 