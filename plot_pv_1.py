""" 
This is an example of how to produce a pv cut from a cube. This follows after the plot_beginner_1.py file.
The script does not rely on any code from etools and thus can be used as is.
There are multiple ways, but the fastest way I found was to use map_coordinates from scipy.ndimage.
The current way implemented in etools uses the RectBivariateSpline from scipy.interpolate but that iterates through each channel map which can be slow at times. 
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import map_coordinates
from scipy.interpolate import interp1d
rad = np.pi / 180. # conversion factor from degree to radians

# set this to the fits file of an image cube
fname = '/lustre/cv/projects/edisk/IRAS04302/export_SBLB/IRAS04302_SBLB_13CO_robust_0.5.pbcor.fits'

# ==== read the data ====
hdul = fits.open(fname)
hdr = hdul[0].header
data = hdul[0].data.T
hdul.close()

# determine the axes
x = (np.arange(hdr['naxis1']) - (hdr['crpix1']-1)) * hdr['cdelt1'] * 3600.
y = (np.arange(hdr['naxis2']) - (hdr['crpix2']-1)) * hdr['cdelt2'] * 3600.

# frequency in Hz
f = (np.arange(hdr['naxis3']) - (hdr['crpix3'] - 1)) * hdr['cdelt3'] + hdr['crval3']

# calculate the velocity axis in km/s
cc = 2.99792458e10
v = cc * (hdr['restfrq'] - f) / hdr['restfrq'] / 1e5

# ==== determine the spatial cut ====
# first, we need to determine the spatial points in the plane of sky

# this is the position axis
p = np.linspace(-5, 5, 1000)

# now determine the location in the sky
# there are multiple ways to determine this. One simple way is to rotate the line by some angle
# note that as implemented here, theta is the angle from the +x direction
# By how the coordinates are setup, the angle is north-of-east, and it's related to the typical position angle (east-of-north) by pi/2 - theta
theta = 45 * rad
xnew = p * np.cos(theta)
ynew = p * np.sin(theta)

# ==== check against the sky image ====
# let's take a look at where we are cutting in the sky plane first

# it's easier to compare against the peak intensity map, but it'll take a few seconds
peak2d = np.max(data, axis=-1)

# now plot
ax = plt.gca()
ax.contourf(x, y, peak2d.T, vmin=0, vmax=0.5 * np.nanmax(peak2d))
ax.plot(xnew, ynew)
plt.show()

# ==== prepare the inputs for map_coordinates ====
# The map_coordinates function treats the inputs in pixel units and not physical units, so we need to convert our physical coordinates to pixel coordinates
xpix = interp1d(x, np.arange(len(x)))(xnew)
ypix = interp1d(y, np.arange(len(y)))(ynew)
zpix = np.arange(len(f), dtype=int)

# calculate the coordinate grids which is required for map_coordinates
# there might be a better way than this
xi, dum = np.meshgrid(xpix, zpix, indexing='ij')
yi, dum = np.meshgrid(ypix, zpix, indexing='ij')
dum, zi = np.meshgrid(np.ones(len(xpix)), zpix, indexing='ij')

# finally we can interpolate
pv = map_coordinates(data, [xi, yi, zi], order=1)

# ==== plotting the pv map ====
ax = plt.gca()
ax.contourf(p, v, pv.T, cmap=plt.cm.RdBu_r, vmin=-np.nanmax(abs(pv)), vmax=np.nanmax(abs(pv)))
plt.show()

