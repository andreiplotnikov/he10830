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
import skimage.transform
import scipy.stats
import scipy.ndimage
import drms
from sunpy.net import Fido
from reproject import reproject_interp
import sunpy.map


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
    aia = a.Instrument('AIA') & a.Wavelength(33 * u.nm, 34 * u.nm)
    
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
    
def limb_darkening(m):
    hpc_coords = sunpy.map.maputils.all_coordinates_from_map(m)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / m.rsun_obs
    r = r.value
    s = m.data
    s = s/np.max(s)
    r = r[s > 0.6]
    s = s[s > 0.6]
    p0 = [0.6, 1, 1, 1, 1]
    plt.plot(r, s, ',')
    plt.plot(r, limb_poly(r, p0))
    #aaa = scipy.optimize.curve_fit(limb_poly, r[r < 1.01], s, p0 = p0)
    k = scipy.optimize.basinhopping(lambda x: np.sum(np.power(limb_poly(r, x) - s, 2)), x0 = p0)
    print(k['fun'])
    aaa = k['x']
    plt.plot(r, limb_poly(r, aaa), ',')
    plt.plot(r, s, ',')
    return(aaa)

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

    fig, ax = plt.subplots()
    
def s0(m):
    hpc_coords = sunpy.map.maputils.all_coordinates_from_map(m)
    r = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) / m.rsun_obs
    r = r.value
    
    return np.median(m.data[r < 0.01])
    
    

directory = 'C:\\data\\bst2\\He10830\\20200923\\'
#directory = 'C:\\data\\bst2\\2\\'

files_list = os.listdir(directory)

# def x2(x, a, k, x0):
#     return a - k*np.power(x - x0, 2)

nodes = 10

file = fits.open(directory + files_list[1])


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

sobel_image = np.hypot(scipy.ndimage.sobel(image, axis=0, mode='constant'),
                       scipy.ndimage.sobel(image, axis=1, mode='constant'))



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
    level = 0.51
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
        print(i)
        
    
    
#plt.pcolormesh(fixed_image)
x = np.array(x)

center_y = image.shape[1]/2
v_line = fixed_image.T[int(center)]
aaa = scipy.interpolate.make_interp_spline(np.arange(len(v_line)), v_line - level*(np.max(v_line) - scatter_level))
means, lengths = pair_means(scipy.interpolate.sproot(aaa))
diameter_y = lengths[np.argmin( np.abs(means - center_y))]
center_y = means[np.argmin( np.abs(means - center_y))]

#fixed_image = np.max(fixed_image)*skimage.transform.rescale(fixed_image/np.max(fixed_image), [diameter/diameter_y, 1])

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
ab_angle = np.arctan(file[0].header.get('ANGLE AB'))



# hdul = fits.HDUList([fits.PrimaryHDU(fixed_image)])
# hdul.writeto('image.fits', overwrite=True)

time = file[0].header.get('TRECSTAR')
timef = time[:4] + '-' + time[5:7]  + '-' + time[8:10] + 'T' + time[11:13] + ':' + time[14:16] + ':' + time[17:19]
time = astropy.time.Time(timef)

coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=time,
                 observer='earth', frame=frames.Helioprojective)

header = sunpy.map.header_helper.make_fitswcs_header(fixed_image, coord,
                                                     reference_pixel=[center, center_y]*u.pixel,
                                                     scale=[scale, scale]*u.arcsec/u.pixel,
                                                     telescope='STT2', rotation_angle = -ab_angle*u.rad - sunpy.coordinates.sun.P(time),
                                                     wavelength = 10830*u.angstrom)

header.update({'STEPSPIX' : file[0].header.get('STEPSPIX'),
               'TRECSTAR' : file[0].header.get('TRECSTAR'),
               'TRECSTOP' : file[0].header.get('TRECSTOP'),
               'OBSERVER' : file[0].header.get('OBSERVER'),
               'ANGLE_AB' : file[0].header.get('ANGLE AB'),
               'SLITWIDT' : file[0].header.get('SLITWIDT'),
               'SLITHEIG' : file[0].header.get('SLITHEIG')})

fixed_image = fixed_image.astype('int16')

fixed_map = sunpy.map.Map(fixed_image, header)

cent_intens = s0(fixed_map)

sunpy.io.fits.write('image.fits', fixed_image, header, overwrite = True)

cor_image = darkening_correction(fixed_map)

cor_map = sunpy.map.Map(cor_image, header)

fixed_map = fixed_map.rotate()

cor_map = cor_map.rotate()

sunpy.io.fits.write('cor_image.fits', cor_image, header, overwrite = True)



fig= plt.figure()

rangex = SkyCoord(-1100 * u.arcsec, -1100 * u.arcsec)
rangey = SkyCoord(1100 * u.arcsec, 1100 * u.arcsec)
fixed_crop = fixed_map.submap(rangex, rangey)
cor_crop = cor_map.submap(rangex, rangey)


ax1 = fig.add_subplot(1, 2, 1, projection=fixed_crop)
fixed_crop.plot(axes = ax1, vmin = 0.5*cent_intens, vmax = 1.1*cent_intens, title = 'Level 1 image')
fixed_crop.draw_limb(axes = ax1)
fixed_crop.draw_grid(axes = ax1)


ax2 = fig.add_subplot(1, 2, 2, projection=cor_crop)
cor_crop.plot(axes = ax2, vmin = 0.8*cent_intens, vmax = 1.1*cent_intens, title = 'Limb-darkening corrected')
cor_crop.draw_limb(axes = ax2)
cor_crop.draw_grid(axes = ax2, grid_spacing = 10 * u.deg)
fig.suptitle(fixed_map.latex_name + ' UTC'  , fontsize = 20, y = 0.9)


plt.savefig('Latest.jpg', bbox_inches = 'tight')
plt.show()



#combined_aia(fixed_map, time)


