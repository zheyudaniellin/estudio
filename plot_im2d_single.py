""" plot_im2d_single.py

Focus on just plotting a single moment image with the continuum. 
Start with the "setting" section under main() and modify the variables for your source. Once that's completed, run this script from command line
python plot_im2d_single.py

"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pdb
import etools as eto
au = eto.natconst.au
rad = np.pi / 180.

def plot_moment(im, cim, mom=0, xlim=None, ylim=None,
    scale_map=1., **kwargs):
    """ plot the moment image
    Parameters
    ----------
    im : eto.image.imageManager
        The image with the moment
    cim : eto.image.imageManager
        the continuum
    mom : int
    xlim, ylim : list of 2 floats
        The x and y limits of the plot in arcseconds, like xlim=[-2,2] and ylim=[-2,2]
    scale_map : float
        Use this to normalize the moment image
    **kwargs are sent to pcolormesh
    """
    ax = plt.gca()

    # plot the colormap
    quant = 'mom%d'%mom
    pc = ax.pcolormesh(im.grid.x, im.grid.y, getattr(im, quant)[:,:,0].T/scale_map, shading='auto', **kwargs)
    cb = plt.colorbar(pc, ax=ax)

    # plot the continuum
    clevs = np.array([5, 10, 50, 100]) * cim.rms
    cc = ax.contour(cim.grid.x, cim.grid.y, cim.I[:,:,0].T, clevs, colors='w', alpha=0.5)

    ax.set_aspect('equal')

    if xlim is not None:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)

    # set up the x-axis formatter to make it show the east direction
    xformatter = FuncFormatter(lambda xval, tpos: '%.1f'%(-xval))
    ax.xaxis.set_major_formatter(xformatter)

    ax.set_xlabel(r'$\Delta$ RA ["]')
    ax.set_ylabel(r'$\Delta$ DEC ["]')

    plt.show()

def main():
    # ==== settings ====
    datdir = 'moments/'
    baseline = 'SB'
    line = '12CO'
    robust = 0.5
    imcor = 'image'
    moment = 8

    # ==== reading ====
    # == moment image ==
    im = eto.image.imageManager()
    fname = eto.kits.get_moment_fits_name(line, baseline, robust, imcor, moment)

    # the moment quantities will be called like 'mom0', 'mom8'
    quant = 'mom%d'%moment
    im.read_fits(os.path.join(datdir, fname), eto.kits.dpc, quant=quant)

    # determine the noise
    # named like 'mom0_rms', 'mom8_rms'
    setattr(im, '%s_rms'%quant, eto.image.get_stats(im, quant, [6, 9])['rms'] )

    # == continuum ==
    cim = eto.kits.easy_read_continuum('SBLB', 0.5)

    # ==== plotting ====
    plot_moment(im, cim, mom=moment, vmin=0)

if __name__ == '__main__':
    main()

