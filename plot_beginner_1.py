""" plot_beginner_1.py
this is an example of how to plot an image cube starting with a fits cube. 
This is independent of other code. 

To use this script, go into python (version 3) and copy and past the commands to learn how to plot things interactively

see also this site for a tutorial
https://learn.astropy.org/tutorials/FITS-images.html
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# set this to the fits file of an image cube
fname = '/lustre/cv/projects/edisk/IRAS04302/export_SBLB/IRAS04302_SBLB_13CO_robust_0.5.pbcor.fits'

# now we want to open the file using the astropy package
hdul = fits.open(fname)

# this is the header
hdr = hdul[0].header

# this is the data
data = hdul[0].data

# usually I like the data in dimensions of x by y by frequency
# so I need to transpose the data
data = data.T

# note that sometimes, the data can have 4 dimensions: x by y by f by Stokes
# since the eDisk project does not have polarization, we usually will only have 2 (continuum or moment maps) or 3 (image cubes) dimensions

# go ahead and close the fits data
hdul.close()

# ======== determining the axis ====
# Let's now determine the x, y, and frequency axes relative to the center
# These are based on the info from the fits header
# There are two places one needs to be careful about:
# - The spatial step size is hdr['cdelt1'] and it is in degrees, but usually we want it in arcseconds. Thus, we need to multiply it by 3600
# - python indexing starts from 0, while fits files start from 1. That's why we need to subtract the hdr['crpix1'] by 1
x = (np.arange(hdr['naxis1']) - (hdr['crpix1']-1)) * hdr['cdelt1'] * 3600.
y = (np.arange(hdr['naxis2']) - (hdr['crpix2']-1)) * hdr['cdelt2'] * 3600.

# Let's get the frequency grid in GHz
f = (np.arange(hdr['naxis3']) - (hdr['crpix3'] - 1)) * hdr['cdelt3'] + hdr['crval3']

# It's usually convenient to look at the velocity grid which requires the rest frequency and the speed of light
# usually the rest frequency is kept in the fits header with hdr['restfreq']
# but sometimes, one may have to define it manually
# speed of light in cm/s
cc = 2.99792458e10

# the velocity here is in km/s
v = cc * (hdr['restfrq'] - f) / hdr['restfrq'] / 1e5

# ==== now try to plot the 2d data ====
""" There are a lot of ways to plot the data, but usually people use the functions imshow, pcolormesh, or contourf from the matplotlib package. 
imshow is often nice, but it can't accept axis arguments which means it only plots in pixel units. Extra work is needed to convert the pixel unit labels to, for example, relative arcseconds. In contrast, pcolormesh and contourf can accept axis arguments. See estudio for examples on how to use imshow and display the correct axis. 

If we want to plot multiple channel maps, then we need to iterate over many indices in the frequency coordinate and plot each slice in different subplots
"""
# get the current plotting axis
ax = plt.gca()

# the data is usually a 3D image cube. To plot a channel map, we need to determine which slice of the cube we want in frequency. This is easily done by indexing in python
# change the last integer into anything you want to view the different channels
pc = ax.pcolormesh(x, y, data[:,:,0].T)

# now we usually want to show it with a colorbar which is easily done by
cb = plt.colorbar(pc, ax=ax)

# one last common procedure is to make sure the aspect ratio is correct (and not deformed)
ax.set_aspect('equal')

plt.show()

# ==== plotting the spectra ====
""" Once we know that the data is a 3D cube, then we know that a spectra is simply the intensity along the frequency axis at some location x, y
"""
# pick whatever index you want for the x and y coordinates
# you will have to think about which index in the x,y coordinates will likely have emission
# Since the disk is usually centered at the image, it is usually at half of the number of spatial pixels of the image
# type data.shape in the command line to see what the dimensions of the data is
plt.plot(f, data[50, 50,:])
plt.show()


# you can also choose to plot the emission as a function of velocity
plt.plot(v, data[50,50,:])
plt.show()

