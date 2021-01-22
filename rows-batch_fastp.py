
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

from tqdm import tqdm

def limb_poly(x, ac):
    a = ac/np.sum(ac)
    s = np.zeros_like(x)
    mu = np.cos(x * np.pi * 0.5)
    for i in range(len(a)):
        s = s + a[i]*np.power(mu, i)
    return s

def pair_means(arr):
    means = []
    lengths = []
    for i in range(len(arr)):
        for j in range(i):
            means.append( 0.5*(arr[i] + arr[j]))
            lengths.append(np.abs(arr[i] - arr[j]))
    return np.array(means), np.array(lengths)

def pearson_shift(image, target, step = 0.1, rng = [-30, 30]):
    table = []
    for shift in np.arange(rng[0], rng[1], step):
        table.append(scipy.stats.pearsonr(scipy.ndimage.shift(image, shift, order = 1),
                                          target)[0])
    return np.arange(rng[0], rng[1], step)[np.argmax(table)]


def limb_poly(x, ac):
    a = ac/np.sum(ac)
    s = np.zeros_like(x)
    mu = np.cos(x * np.pi * 0.5)
    for i in range(len(a)):
        s = s + a[i]*np.power(mu, i)
    return s
            

def robust_std(a):
    arr = a.copy()
    arr -= np.median(arr)
    arr = np.abs(shift_list)
    return np.quantile(arr, 0.667)

def combined_aia(fixed_map, time):
    atime = a.Sample(24 * u.hour) & a.Time(time, time + astropy.time.TimeDelta(1e-3)) 
    aia = a.Instrument.aia & a.Wavelength(17 * u.nm, 18 * u.nm)
    
    res = Fido.search(atime, aia)
    files = Fido.fetch(res[:, 0])
    
    map_aia = sunpy.map.Map(files).resample((1024, 1024)*u.pix)
    output, footprint = reproject_interp(fixed_map, map_aia.wcs, map_aia.data.shape)
    
    out_fixed = sunpy.map.Map(output, map_aia.wcs)
    
    fig = plt.figure(figsize=(20, 15))
    ax1 = fig.add_subplot(1, 1, 1, projection=map_aia)
    map_aia.plot(axes=ax1)
    out_fixed.plot(axes=ax1, alpha=0.5)
    
    plt.show()
    
    plt.plot(map_aia.data.flatten(), out_fixed.data.flatten(), ',')
    
def circle_check(im, x0, y0, d, i, step = 0.1, rng = [-5, 5]):
    image = im.copy()
    
    y_c, x_c = np.indices(image.shape).astype('float64')
    x_c -= x0
    y_c -= y0
    r = np.sqrt(x_c ** 2 + y_c ** 2) / (0.5*d)
    
    table = []
    
    for shift in np.arange(rng[0], rng[1], step):
        image[i] = scipy.ndimage.shift(im[i], shift, order = 1)
        table.append(np.sum(image[r < 1])/np.sum(im))
        
    return np.arange(rng[0], rng[1], step)[np.argmax(table)]

def limb_darkening(m):
    hpc_coords = sunpy.map.maputils.all_coordinates_from_map(m)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / m.rsun_obs
    plt.plot(r.flatten(), m.data.flatten(), ',')
    
def heliographic(m):
    shape_out = [720, 1440]
    frame_out = SkyCoord(0, 0, unit=u.deg,
                         frame="heliographic_stonyhurst",
                         obstime=m.date)
    header = sunpy.map.make_fitswcs_header(shape_out,
                                           frame_out,
                                           scale=[180 / shape_out[0],
                                                  360 / shape_out[1]] * u.deg / u.pix,
                                           projection_code="CAR")
    
    out_wcs = astropy.wcs.WCS(header)
    array, footprint = reproject_interp(m, out_wcs, shape_out=shape_out)
    outmap = sunpy.map.Map((array, header))
    outmap.plot_settings = m.plot_settings    
    
    fig = plt.figure()
    ax = plt.subplot(projection=outmap)
    outmap.plot(ax)
    
    ax.set_xlim(0, shape_out[1])
    ax.set_ylim(0, shape_out[0])
    
    plt.show()
    
def darkening_correction(m):
    
    coeffs = [ 0.64434808,  1.29852297, -2.62043326,  2.72090958, -1.04334737]
    
    hpc_coords = sunpy.map.maputils.all_coordinates_from_map(m)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / m.rsun_obs
    r = r.value
    s = m.data
    mask = np.where(r < 1, limb_poly(r, coeffs), 1)
    s = s/mask
    
    s = s.astype('int16')
    
    return s

def s0(m):
    hpc_coords = sunpy.map.maputils.all_coordinates_from_map(m)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / m.rsun_obs
    r = r.value
    
    return np.median(m.data[r < 0.01])
    

def time_from_file(file_name):
    file = fits.open(file_name)
    time = file[0].header.get('TRECSTAR')
    year = int(time[:4])
    year_short = time[2:4]
    month = int(time[5:7])
    day = int(time[8:10])
    hour = int(time[11:13])
    minute = int(time[14:16])
    sec = int(time[17:19])
    
    timezone = astropy.time.TimezoneInfo(utc_offset=0*u.hour)
    dt = datetime.datetime(year, month, day, hour, minute, sec, tzinfo=timezone)
    
    time = astropy.time.Time(dt)   
    
    return time

    

directory = 'C:\\data\\bst2\\He10830\\'

#directory = 'D:\\20200731\\'


save_files = False

save_path = 'C:\\results\\helium\\'

files_list = []

for root, dirs, files in os.walk(directory, topdown = False):
    for name in files:
        files_list.append(os.path.join(root, name))
        
obs_except = 0
proceeded = 0

start_time = time_from_file(files_list[0])
finish_time = time_from_file(files_list[-2])

nodes = 100

time_linspace = astropy.time.Time(list(start_time + i*(finish_time - start_time)/nodes for i in range(nodes)))
delta_linspace = list( (i*(finish_time - start_time)/nodes).value for i in range(nodes))

angles = sunpy.coordinates.sun.P(time_linspace)

solar_p_spline = scipy.interpolate.make_interp_spline(delta_linspace, angles)


for file_name in tqdm(files_list):
    try:
        time = time_from_file(file_name)
        
        file = fits.open(file_name)
        
        time = file[0].header.get('TRECSTAR')
        year = int(time[:4])
        year_short = time[2:4]
        month = int(time[5:7])
        day = int(time[8:10])
        hour = int(time[11:13])
        minute = int(time[14:16])
        sec = int(time[17:19])
        
        if save_files:
            try:
                os.mkdir(save_path + '%s_%02i' % (year_short, month))
            except:
                pass
            local_save = save_path + '%s_%02i\\' % (year_short, month)
            
        image = file[0].data
        
        center = image.shape[0]/2
        
        fixed_image = image.copy()
        lines = np.empty_like(image)
        
        scale = 2.97*701/image.shape[0]
        
        x = []
        
        shift_list = []
        
        scatter_level = np.min(image[image > 0])
        
        #image = skimage.filters.gaussian(image, 3)
        
        diameter = 0
    
        for i in range(image.shape[0]):
        #for i in range(30, 180):
            line = image[i]
            # line = line[200:500]
            # vec = scipy.optimize.curve_fit(x2, np.arange(200, 500), line)[0]
            # x0 = vec[2]
            # x.append(vec)
            # shift = 0.5*(pearson_shift(image[i], fixed_image[i - 1]) + 
            #             pearson_shift(image[i], fixed_image[i - 1]))
            # lines[i] = x2(np.arange(0, image.shape[1]), vec[0], vec[1], vec[2] - shift)
            level = 0.5
            if np.max(image[i]) > 2*scatter_level:
                local_scatter = np.min(line[line > 0])
                aaa = scipy.interpolate.make_interp_spline(np.arange(image.shape[1]), line - level*(np.max(line) - scatter_level))       
                means, lengths = pair_means(scipy.interpolate.sproot(aaa))
                shift = means[np.argmin( np.abs(means - center))] - center
    
                    
                diameter = np.max( (lengths[np.argmin( np.abs(means - center))], diameter))
            
                shift_list.append(shift)
            
                fixed_image[i] = scipy.ndimage.shift(image[i], -shift, order = 1)
            else:
                fixed_image[i] = image[i]
                shift_list.append(0)
            #print(i)
        
        rstd = robust_std(shift_list)
        shift_med = np.median(shift_list)
        
        # plt.pcolormesh(fixed_image)
        # plt.show()
        
        for i in range(image.shape[0]):
            if np.abs(shift_list[i] - shift_med) > 5*rstd:
                shift = 0.5*(pearson_shift(image[i], fixed_image[i - 1]) + 
                          pearson_shift(image[i], fixed_image[i + 1]))
                fixed_image[i] = scipy.ndimage.shift(image[i], -shift, order = 1)
                #print(i)
            
        
        
        #plt.pcolormesh(fixed_image)
        x = np.array(x)
        
        center_y = image.shape[1]/2
        v_line = fixed_image.T[int(center)]
        aaa = scipy.interpolate.make_interp_spline(np.arange(len(v_line)), v_line - level*(np.max(v_line) - scatter_level))
        means, lengths = pair_means(scipy.interpolate.sproot(aaa))
        diameter_y = lengths[np.argmin( np.abs(means - center_y))]
        center_y = means[np.argmin( np.abs(means - center_y))]
        
        if diameter/diameter_y > 0.5 and diameter/diameter_y < 1.5:
            fixed_image = np.max(fixed_image)*skimage.transform.rescale(fixed_image/np.max(fixed_image), [diameter/diameter_y, 1])
        
        v_line = fixed_image.T[int(center)]
        aaa = scipy.interpolate.make_interp_spline(np.arange(len(v_line)), v_line - level*(np.max(v_line) - scatter_level))
        means, lengths = pair_means(scipy.interpolate.sproot(aaa))
        diameter_y = lengths[np.argmin( np.abs(means - center_y))]
        center_y = means[np.argmin( np.abs(means - center_y))]
        
        # theta = np.linspace(0, 2*np.pi, 100)
        # x1 = diameter/2*np.cos(theta) + center
        # x2 = diameter/2*np.sin(theta) + center_y
        
        # plt.pcolormesh(fixed_image)
        # plt.plot(x1, x2, linewidth = 5)
        # plt.show()
        
        
        fixed_image = np.flip(fixed_image, axis = 1)
        ab_angle = -np.arctan(file[0].header.get('ANGLE AB'))
    
        # hdul = fits.HDUList([fits.PrimaryHDU(fixed_image)])
        # hdul.writeto('image.fits', overwrite=True)
    
    
    
        coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=time,
                         observer='earth', frame=frames.Helioprojective)
        
        p_angle = solar_p_spline( (time - start_time).value)
    
        header = sunpy.map.header_helper.make_fitswcs_header(fixed_image, coord,
                                                             reference_pixel=[center, center_y]*u.pixel,
                                                             scale=[scale, scale]*u.arcsec/u.pixel,
                                                             rotation_angle = ab_angle*u.rad - p_angle,
                                                             telescope='STT2', observatory = 'CrAO', 
                                                             wavelength = 10830*u.angstrom)
        
        header.update({'STEPSPIX' : file[0].header.get('STEPSPIX'),
               'TRECSTAR' : file[0].header.get('TRECSTAR'),
               'TRECSTOP' : file[0].header.get('TRECSTOP'),
               'ANGLE_AB' : file[0].header.get('ANGLE AB'),
               'SLITWIDT' : file[0].header.get('SLITWIDT'),
               'SLITHEIG' : file[0].header.get('SLITHEIG')})
    
        try:
            header.update({'OBSERVER' : file[0].header.get('OBSERVER')})
        except:
            obs_except += 1
            pass
            
        
        fixed_image = fixed_image.astype('int16')
        
        fixed_map = sunpy.map.Map(fixed_image, header)
        
        cent_intens = s0(fixed_map)
        
        cor_image = darkening_correction(fixed_map)       
        cor_map = sunpy.map.Map(cor_image, header)
        
        fixed_map = fixed_map.rotate()
        cor_map = cor_map.rotate()
        
        rangex = SkyCoord(-1100 * u.arcsec, -1100 * u.arcsec)
        rangey = SkyCoord(1100 * u.arcsec, 1100 * u.arcsec)
        
        fixed_map = fixed_map.submap(rangex, rangey)
        cor_map = cor_map.submap(rangex, rangey)
        

    except:
        continue
    
    if save_files:
        sunpy.io.fits.write(local_save + '%s%02i%02i_%02i%02i.fits' % (year_short, month, day, hour, minute), fixed_map.data.astype('int16'), fixed_map.fits_header, overwrite = True)
    
        try:
            fig= plt.figure()
            
            ax1 = fig.add_subplot(1, 2, 1, projection=fixed_map)
            fixed_map.plot(axes = ax1, vmin = 0.5*cent_intens, vmax = 1.1*cent_intens)
            fixed_map.draw_limb(axes = ax1)
            fixed_map.draw_grid(axes = ax1, grid_spacing = 10 * u.deg)
            
            
            ax2 = fig.add_subplot(1, 2, 2, projection=cor_map)
            cor_map.plot(axes = ax2, vmin = 0.8*cent_intens, vmax = 1.1*cent_intens)
            cor_map.draw_limb(axes = ax2)
            cor_map.draw_grid(axes = ax2, grid_spacing = 10 * u.deg)
            
            fig.suptitle(fixed_map.latex_name + ' UTC'  , fontsize = 20, y = 0.9)
            
            plt.savefig(local_save + '%s%02i%02i_%02i%02i.jpg' % (year_short, month, day, hour, minute), bbox_inches = 'tight')
            plt.close()
            proceeded += 1
        except:
            pass

#combined_aia(fixed_map, time)


