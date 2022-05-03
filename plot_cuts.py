""" plot_cuts.py
produce cuts along the continuum image
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
import astropy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pdb
import etools as eto
au = eto.natconst.au
rad = np.pi / 180.
figlab = 'abcdefghijklmnpqrstuvwxyz'
lincol = ['#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f']

def plot_sky(im, cuts, tags):
    """ visualize where the cuts are taken in the sky compared to the image
    """
    fig = plt.figure()
    ax = fig.gca()
    pc = ax.contourf(im.grid.x, im.grid.y, im.I[:,:,0].T)

    # put some contours too
    clevs = np.array([3]) * im.rms
    cc = ax.contour(im.grid.x, im.grid.y, im.I[:,:,0].T, clevs, colors='w')

    ax.set_aspect('equal')

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Plot beam
    xc = xlim[0] + 0.80 * np.diff(xlim)
    yc = ylim[0] + 0.10 * np.diff(ylim)
    im.plot_beam(ax=ax, axis_unit='arcsec', beamxy=(xc, yc))

    # Plot the cuts
    for i, icut in enumerate(cuts):
        ax.plot(icut.xpnt, icut.ypnt, label=tags[i], color='C%d'%i)

        # also plot the arrow
        dx = np.diff(icut.xpnt)
        dy = np.diff(icut.ypnt)
        ax.arrow(icut.xpnt[-1], icut.ypnt[-1], dx[-1], dy[-1],
            shape='full', lw=10, color='C%d'%i,
            length_includes_head=True, head_width=.05)

    ax.legend()

    fig.tight_layout()
    plt.show()

def plot_simple(cuts, tags, figsize=None, pdfname=None):
    """ plot the cuts in mjy/beam
    cuts : list of eto.image.imageCut objects
    tags : list of str
        label for each cut
    figsize : tuple of 2 floats
        the size of the figure
    pdfname : str
        the name of the pdf file
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.gca()

    for i, icut in enumerate(cuts):
        ax.plot(icut.laxis, icut.I[:,0]/1e-3, label=tags[i], linewidth=2, color=lincol[i])

    ax.legend()
    ax.set_xlabel('location ["]')
    ax.set_ylabel(r'$I_{\nu}$ [mJy/beam]')

    fig.tight_layout()
    if pdfname is not None:
        fig.savefig(pdfname)

    plt.show()

def plot_dual_axis(cuts, tags, thres=3, tb_ticks=None,
    cols=None, do_beam=False,
    pdfname=None, figsize=(10,5)):
    """ plot the cuts in mjy/beam and Tb and also arcsec and au
    all just in one plot
    The cuts must have the same beam 
    Use this with caution. It's not quite well tested yet. 
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.gca()

    # basic stuff
    bmaj = cuts[0].psfarg['bmaj']
    bmin = cuts[0].psfarg['bmin']
    wav = cuts[0].w

    # plot all the cuts first
    if cols is None:
        cols = lincol

    xmax = 0
    for i, icut in enumerate(cuts):
        ax.plot(icut.laxis, icut.I[:,0]/1e-3, label=tags[i], linewidth=2, color=cols[i])
        if icut.laxis.max() > xmax:
            xmax = icut.laxis.max()

    # ==== set the xlim ====
    # use the most extended among the cuts
    ax.set_xlim(-xmax, xmax)

    # ==== calculate the au labels ====
    dpc = cuts[0].dpc
    def au_to_arcsec(x):
        return x / dpc
    def arcsec_to_au(x):
        return x * dpc
    ax_top = ax.secondary_xaxis('top', functions=(arcsec_to_au, au_to_arcsec))

    ax.set_xlabel('location ["]')
    ax_top.set_xlabel('location [au]')

    # ==== calculate the Tb labels ====
    # it's easier to set the ylimit down to 0
#    ax.set_ylim(0, None)

    # prepare the conversion between jy/beam to brightness temperature
    def mjypbeam_to_tb(x):
        cgs = eto.image.jypbeam_to_cgs(x/1e3, bmaj, bmin)
        tb = eto.image.cgs_to_tb(cgs, wav)
        return tb

    def tb_to_mjypbeam(x):
        cgs = eto.image.tb_to_cgs(x, wav)
        jypbeam = eto.image.cgs_to_jypbeam(cgs, bmaj, bmin) / 1e-3
        return jypbeam

    ax_right = ax.secondary_yaxis('right', functions=(mjypbeam_to_tb, tb_to_mjypbeam))

    ax.set_ylabel(r'$I_{\nu}$ [mjypbeam]')
    ax_right.set_ylabel(r'$T_{b}$ [K]')

    # ==== beam ====
    if do_beam:
        fwhm = np.sqrt(bmaj*bmin)
        beamx = np.linspace(-2*fwhm, 2*fwhm, 129)
        beamy = fn_gaussian(beamx, cuts[0].I.max()/1e-3, 0, fwhm)
        ax.plot(beamx, beamy, color='k', linestyle=':', label='Beam')

    # ==== shade the noise region ====
    rms = cuts[0].rms
    xlim = ax.get_xlim()
    ax.fill_between(xlim, [0,0], [thres*rms/1e-3]*2, color='grey', alpha=0.3)

    ax.legend()

    fig.tight_layout()
    if pdfname is not None:
        fig.savefig(pdfname)

    plt.show()

def main():
    # ==== settings ====
    baseline = 'SBLB'
    robust = 0.5

    # the disk position angle of the major axis (east of north)
    diskpa = -5 * rad

    # continuum
    cim = eto.kits.easy_read_continuum(baseline, robust)

    # the wavelength is necessary because it's not given for continuum
    cim.grid.set_wavelength(np.array([1300]))
    cim.grid.get_frequency()

    # trim the image, so that it's easier for interpolation
    trimpar = {'xlim':[-2.5, 2.5], 'ylim':[-2.5, 2.5]}
    cim.trim(trimpar)

    # make step size
    stepsize = np.sqrt(cim.psfarg['bmaj']*cim.psfarg['bmin']) / 8.

    # ==== get the major axis ====
    # laxis is the path array 
    laxis = np.arange(-2.3, 2.3, stepsize)
    trackkw = {'x0':0, 'y0':0, 'theta':(diskpa)}
    major = cim.cut(laxis, quant=['I'], track='linear', trackkw=trackkw)

    # get the minor axis
    laxis = np.arange(-0.45, 0.4, stepsize)
    trackkw['theta'] = trackkw['theta'] + 90 * rad
    minor = cim.cut(laxis, trackkw=trackkw)

    # keep some of the properties from the image for the cuts
    for ikey in ['rms', 'psfarg']:
        setattr(major, ikey, getattr(cim, ikey))
        setattr(minor, ikey, getattr(cim, ikey))

    # ==== plotting ====
    # plot where the cuts are taken in the sky
    plot_sky(cim, [major, minor], ['major', 'minor'])

    # plot the profiles of the cuts
    plot_simple([major, minor], ['major', 'minor'])

    # use this with caution. It's not quite well tested yet
#    pdfname = 'results/cont_cuts.pdf'
    plot_dual_axis([major, minor], ['major', 'minor'], do_beam=False, figsize=(12,5), pdfname=None)


if __name__ == '__main__':
    main()

