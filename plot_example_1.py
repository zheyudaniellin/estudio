""" plot_example_1.py
This is the simplest way to plot something interactively and has the least amount of dependencies. 
Change the parameters, 'dpc' and 'fname', to match your data.
Then, simply go into python 3 from your terminal and then start copy and pasting the commands
"""
import numpy as np
import matplotlib.pyplot as plt
import etools as eto

# distance to source in pc
dpc = 140.

# determine the fits name of the image you want to plot
# it can be the continuum or the moment images or an image cube
fname = '/lustre/cv/projects/edisk/IRAS04302/export_SBLB/IRAS04302_SBLB_continuum_robust_0.0.pbcor.tt0.fits'

# prepare the imageManager
im = eto.image.imageManager()

# read the fits image
im.read_fits(fname, dpc, quant='I')

"""
Notice that the imageManager now has an attribute 'I' which is a 3d numpy array in x by y by frequency and an attribute 'grid' which is the gridManager and contains info about the grid
"""

# start plotting 
ax = plt.gca()

pc = ax.pcolormesh(im.grid.x, im.grid.y, im.I[:,:,0].T, cmap=plt.cm.gist_heat)

cb = plt.colorbar(pc, ax=ax)

ax.set_aspect('equal')

plt.show()

"""
1. For 2D images, like the continuum or moment maps, the third index should be 0, like im.I[:,:,0]
But if you have a channel cube, then use the third index to choose your frequency, like im.I[:,:,50] 

The frequency grid is kept as: 
im.grid.f

If you want to get the velocity grid from the frequency, simply do: 
im.grid.get_velocity()

Note that this requires you to set the restfrequency, which is usually done for you if the fits image is a channel cube. If not, then you'll want to set this before using get_velocity: 
im.grid.restfreq = the rest frequency in Hz

2. Note that the x-axis of the grid is towards WEST. This made life a lot easier for me. If you want to show the x-axis directed towards east, then look at plot_im2d_single.py for example on using an axis formatter. 
"""
