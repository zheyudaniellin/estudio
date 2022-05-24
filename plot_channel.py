""" plot_channel.py
use this script to plot channel maps
"""
import numpy as np
import matplotlib.pyplot as plt
import pdb
import pickle
import os
import etools as eto
au = eto.natconst.au
rad = np.pi / 180.

def plot_basic(im, nrow=5, cim=None, vrange=[-5, 15], perc=[10, 90],
    xlim=None, ylim=None, figsize=None, axis_ticks=None,
    clevs=None,):
    """ plot the channels
    nrow : int
        the number of rows to be plotted
    cim : eto.image.imageManager
        the continuum if desired
    vrange : list
        the range of velocity in km/s
    """
    # basic common info
    pixel_size = abs(im.grid.dx)

    # determine velocity range
    reg = (vrange[0] <= im.grid.v) & (im.grid.v <= vrange[1])
    vel = im.grid.v[reg]
    vinx = np.arange(len(im.grid.v))[reg]
    nchan = len(vel)
    v_a, v_b = vinx[0], vinx[-1]

    # determine the intensity range
    if xlim is not None:
        x_a = np.argmin(abs(im.grid.x-xlim[0]))
        x_b = np.argmin(abs(im.grid.x-xlim[1]))
    else:
        x_a, x_b = 0, im.grid.nx-1
    xpixel_line = np.arange(len(im.grid.x[x_a:x_b]))

    if ylim is not None:
        y_a = np.argmin(abs(im.grid.y - ylim[0]))
        y_b = np.argmin(abs(im.grid.y - ylim[1]))
    else:
        y_a, y_b = 0, im.grid.ny-1
    ypixel_line = np.arange(len(im.grid.y[y_a:y_b]))

#    vlim = np.percentile(im.I[x_a:x_b,y_a:y_b,v_a:v_b].flatten(), perc)
    vlim = [None, None]

    # determine the spatial ticks
    if axis_ticks is None:
        axis_ticks = np.array([-3, -2, -1, 0, 1, 2, 3])

    # prepare the continuum
    # need the x,y axis on the line pixel grid for contour
    if cim is not None:
        xpixel_cont = np.interp(cim.grid.x, im.grid.x[x_a:x_b], xpixel_line)
        ypixel_cont = np.interp(cim.grid.y, im.grid.y[y_a:y_b], ypixel_line)

    # ==== plotting ====
#    nrow = np.floor(np.sqrt(nchan))
    ncol = np.ceil(nchan / nrow)
    nrow, ncol = int(nrow), int(ncol)
    fig, axgrid = plt.subplots(nrow,ncol,sharex=True,sharey=True,
        squeeze=False, figsize=figsize)
    axes = axgrid.flatten()

    # iterate through each image
    for i in range(nchan):
        ax = axes[i]

        # plot the line image
        ax.imshow(im.I[x_a:x_b,y_a:y_b,vinx[i]].T, origin='lower',
            cmap=plt.cm.gist_heat,
            vmin=vlim[0], vmax=vlim[1])

        # plot continuum if desired
        # skip this for now
        if cim is not None:
#            levs = np.array([0.1]) * np.nanmax(cont.I)
            levs = np.array([5]) * cim.rms
            ax.contour(xpixel_cont, ypixel_cont, cim.I[:,:,0].T,
                levs, colors='w',alpha=0.3)

        ax.set_aspect('equal')

        # ==== other annotations ====
        # velocity text
        txt = '%.2f'%vel[i]
        ax.text(0.05, 0.95, txt, transform=ax.transAxes, color='w',
            va='top', ha='left')

        # plot beam
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xc = xlim[0] + 0.1 * np.diff(xlim)
        yc = ylim[0] + 0.1 * np.diff(ylim)
#        im.plot_beam(ax=ax, beamxy=(xc,yc), facecolor='w', axis_unit='arcsec')
        eto.image.plot_beam(ax, xc, yc, im.psfarg['bmaj']/pixel_size, im.psfarg['bmin']/pixel_size, im.psfarg['bpa'])

        # plot the length
        xc = xlim[0] + 0.9 * np.diff(xlim)
        yc = ylim[0] + 0.1 * np.diff(ylim)
        eto.image.plot_length(ax, xc, yc, 100 / im.grid.dpc / pixel_size, '')

        # add the tick labels
        xpixel_cen = len(im.grid.x[x_a:x_b])//2
        ypixel_cen = len(im.grid.y[y_a:y_b])//2
        ax.set_xticks(xpixel_cen + axis_ticks / pixel_size)
        ax.set_yticks(ypixel_cen + axis_ticks / pixel_size)
        ax.tick_params(axis='both', direction='in', color='w')

        ax.set_xticklabels([])
        ax.set_yticklabels([])

        # Make the axes white as well
        for side in ["left","right","top","bottom"]:
            ax.spines[side].set_color("white")

#        ax.set_axis_off()

    if nchan != (nrow*ncol):
        for i in range(nchan, nrow*ncol, 1):
            axes[i].set_axis_off()

#    fig.tight_layout(w_pad=0.05, h_pad=0.05)
    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.show()

def main():
    # ==== settings ====
    # distance to source in pc
    dpc = 140.

    # ==== read channel map ====
    # this is the name of the fits file
    fname = '/lustre/cv/projects/edisk/IRAS04302/export_SBLB/IRAS04302_SBLB_13CO_robust_0.5.image.fits'
    # prepare the image object
    im = eto.image.imageManager()

    # read the fits file and also set the distance
    im.read_fits(fname, dpc)

    # set the rest frequency
    # estudio keeps the header of the fits file as an attribute of grid
    im.grid.set_restfreq(im.grid.fits_hdr['restfrq'])

    # calculate the velocity
    im.grid.get_velocity()

    # if we want to trim the image
    par = {'xlim':[-10,10], 'ylim':[-10, 10]}
    im.trim(par)

    # if we want to average the channels
#    im.average_channel(3)

    # ==== read continuum ====
    cim = eto.kits.easy_read_continuum('SBLB', 0.5)
    cim.trim({'xlim':[-2.5, 2.5], 'ylim':[-2.5,2.5]})

    # ==== plotting ====
    vsys = 5.9
    delv = 6.
    axis_ticks = np.arange(-5, 5.01, 2.5)
    plot_basic(im, cim=cim, nrow=6, axis_ticks=axis_ticks,
        xlim=[-6,6], ylim=[-6, 6], vrange=[vsys-delv, vsys+delv])

if __name__ == '__main__':
    main()


