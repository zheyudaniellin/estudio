""" plot_im2d_multi.py
plot multiple images
The most important function here is plot_moment_cmap
You can change the coloring scheme and other settings there, so be sure to be very familiar with it
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
import pdb
import etools as eto
au = eto.natconst.au
rad = np.pi / 180.

def plot_moment_cmap(cim, ims, lines, mom=0,
    mask_mom=None, mthres=3, mground=0.,
    gridmode='physical',
    fig=None, axes=None, noshow=False,
    figsize=None, xlim=None, ylim=None, vmin=None, vmax=None, **kwargs):
    """ plot multiple moment0 as color maps
    Parameters
    ----------
    cim : eto.image.imageManager
    ims : list of eto.image.imageManager
        the images
    lines : list of str
        names of the lines
    mom : int
        moment to be plotted as the colormap
    mask_mom : int
        if None, then don't consider masking
    mthres : float
        some multiplication of rms
    mground : float
        som multiplication of avg

    gridmode : str
        The mode to treat the grid. This determines if we use pcolormesh or imshow. Using imshow makes it easier to save to pdf files
        'physical' = use arcseconds which allows pcolormesh/contourf
        'pixel' = use pixels which allows imshow
    xlim, ylim : list of 2 floats
        the values in arcseconds (even if gridmode='pixel')

    noshow : bool
        if True, then don't show the plot in the end. This is useful for saviging it into a file
    fig : matplotlib figure
        if not given, then it will be setup from plt.subplots
    axes : list of matplotlib axes
        if not given, then it will be setup from plt.subplots
    figsize : tuple of 2 floats
        size of the figure if setting up from plt.subplots
    vmin, vmax : float
        The value minimum and maximum of the color maps
    other **kwargs will be passed to the colormap plotting routine like pcolormesh/contourf/imshow
    """
    # determine if the figure object and ax for plots are provided
    # if provided, then use those. else create subplots
    if fig is None:
        nrow, ncol = 1, len(ims)
        fig, axgrid = plt.subplots(nrow,ncol,sharex=True,sharey=True,
            squeeze=False, figsize=figsize)
        axes = axgrid.flatten()
    else:
        # check if the number of ax are the same as number of images
        if len(axes) != len(ims):
            raise ValueError('number of ax does not match number of images')

    for i in range(len(ims)):
        ax = axes[i]
        im = ims[i]

        ax.set_facecolor('w')

        # choose the 2d image and scale based on the units
        if mom == 0:
            im2d = im.mom0[:,:,0] * 1e3
        elif mom == 1:
            im2d = im.mom1[:,:,0] * 1.
        elif mom == 8:
            im2d = im.mom8[:,:,0] * 1e3
        elif mom == 9:
            im2d = im.mom9[:,:,0] * 1.
        elif mom == 10:
            im2d = im.mom10[:,:,0] * 1e3
        else:
            raise ValueError('moment unknown')
        # clip the limits
        if (vmin is not None) or (vmax is not None):
            im2d = np.clip(im2d, vmin, vmax)

        # mask it if desired
        if mask_mom is not None:
            mask = getattr(im, 'mom%d'%mask_mom)[:,:,0]
            thres = getattr(im, 'mom%d_rms'%mask_mom) * mthres
            if mground != 0:
                thres += getattr(im, 'mom%d_avg'%mask_mom) * mground
            reg = mask <= thres
            im2d[reg] = np.nan

        # plot the colormap
        if gridmode == 'physical':
            pc = ax.contourf(im.grid.x, im.grid.y, im2d.T, 32,
                vmin=vmin, vmax=vmax, **kwargs)
        elif gridmode == 'pixel':
            pc = ax.imshow(im2d.T, origin='lower', vmin=vmin, vmax=vmax, **kwargs)
        else:
            raise ValueError('gridmode unknown')

        # colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="2%", pad=0.5)
        cb = fig.colorbar(pc, cax=cax, orientation='horizontal')

        # Plot the continuum
        clevs = np.array([5, 25, 50, 100]) * cim.rms
        ckwargs = {'colors':'k', 'alpha':0.5}
        if gridmode == 'physical':
            cc = ax.contour(cim.grid.x, cim.grid.y, cim.I[:,:,0].T, clevs, **ckwargs)
        elif gridmode == 'pixel':
            xpixel_line = np.arange(im.grid.nx)
            ypixel_line = np.arange(im.grid.ny)
            xpixel_cont = np.interp(cim.grid.x, im.grid.x, xpixel_line)
            ypixel_cont = np.interp(cim.grid.y, im.grid.y, ypixel_line)
            cc = ax.contour(xpixel_cont, ypixel_cont, cim.I[:,:,0].T, clevs, **ckwargs)
        else:
            raise ValueError('gridmode unknown')

        ax.set_aspect('equal')
        # set the x and y limits
        if xlim is not None:
            if gridmode == 'physical':
                xlimpar = xlim
            elif gridmode == 'pixel':
                xlimpar = np.interp(xlim, im.grid.x, xpixel_line)
            ax.set_xlim(xlimpar)

        if ylim is not None:
            if gridmode == 'physical':
                ylimpar = ylim
            elif gridmode == 'pixel':
                ylimpar = np.interp(ylim, im.grid.y, ypixel_line)
            ax.set_ylim(ylimpar)

        # plot the length
        xlimpar = ax.get_xlim()
        ylimpar = ax.get_ylim()
        xc = xlimpar[0] + 0.80 * np.diff(xlimpar)
        yc = ylimpar[0] + 0.05 * np.diff(ylimpar)
        if gridmode == 'physical':
            pltlen = 100. / im.grid.dpc
        elif gridmode == 'pixel':
            pltlen = 100. / abs(im.grid.dx * im.grid.dpc)
        eto.image.plot_length(ax, xc, yc, pltlen, '100 au', color='k')

        # plot the beam
        xc = xlimpar[0] + 0.2 * np.diff(xlimpar)
        yc = ylimpar[0] + 0.05 * np.diff(ylimpar)
        if gridmode == 'physical':
            bscale = 1.
        elif gridmode == 'pixel':
            bscale = abs(im.grid.dx)
        eto.image.plot_beam(ax, xc, yc, im.psfarg['bmaj']/bscale, im.psfarg['bmin']/bscale, im.psfarg['bpa'], facecolor='k')

        # we want to plot in physical arcsecond units and not pixel units
        if gridmode == 'pixel':
            xformatter = FuncFormatter(lambda xval, tpos: '%.1f'%(-np.interp(xval, xpixel_line, im.grid.x )) )
            ax.xaxis.set_major_formatter(xformatter)
            yformatter = FuncFormatter(lambda xval, tpos: '%.1f'%np.interp(xval, ypixel_line, im.grid.y ) )
            ax.yaxis.set_major_formatter(yformatter)

        # annotate stuff
        ax.set_title(lines[i])

    fig.tight_layout()
    if not noshow:
        plt.show()

def plot_one_at_a_time(cim, ims, line, pdfname, mom=0, **kwargs):
    """ convenient function to iterate through the lines to
    produce one pdf file of continuum with lines one at a time
    """
    pdf = PdfPages(pdfname)

    for i in range(len(ims)):
        fig = plt.figure()
        plot_moment_cmap(cim, [ims[i]], [line[i]], mom=mom, noshow=True, fig=fig, axes=[fig.gca()], gridmode='pixel', **kwargs)
        pdf.savefig(fig)
        plt.close()
    pdf.close()

def main():
    # ==== settings ====
    src = [
        # line, robust, baseline
        ['12CO', 0.5, 'SB'],
        ['13CO', 0.5, 'SB'],
        ['C18O', 0.5, 'SB'],
        ['SO', 0.5, 'SB'],
#        ['H2CO_3_21-2_20_218.76GHz', 0.5, 'SBLB'],
#        ['H2CO_3_03-2_02_218.22GHz', 0.5, 'SBLB'],
#        ['H2CO_3_22-2_21_218.47GHz', 0.5, 'SBLB'],
#        ['c-C3H2_217.82', 0.5, 'SBLB'],
#        ['cC3H2_217.94', 0.5, 'SBLB'],
#        ['cC3H2_218.16', 0.5, 'SBLB'],
#        ['DCN', 0.5, 'SBLB'],
#        ['CH3OH', 0.5, 'SBLB'],
#        ['SiO', 0.5, 'SBLB'],
        ]
    imcor = 'image'
    trimpar = {'xlim':[-15, 15], 'ylim':[-15,15]}

    # choose the moment to be plotted
    moment = 9

    # ==== orgainze the lines ====
    line = [isrc[0] for isrc in src]
    robust = [isrc[1] for isrc in src]
    baseline = [isrc[2] for isrc in src]

    # ==== reading ====
    # read lines
    ims = []
    for i in range(len(src)):
        im = eto.image.imageManager()
        fname = eto.kits.get_moment_fits_name(line[i], baseline[i], robust[i], imcor, moment)
        fname = os.path.join('moments', fname)
        iquant = 'mom%d'%moment
        im.read_fits(fname, eto.kits.dpc, quant=iquant)
        stats = eto.image.get_stats(im, iquant, [6, 9])
        setattr(im, '%s_rms'%iquant, stats['rms'] )
        setattr(im, '%s_avg'%iquant, stats['avg'] )
#        im.trim(trimpar)
        ims.append(im)

    # read continuum
    cim = eto.kits.easy_read_continuum('SBLB', 0.5)
    cim.trim(trimpar)

    # ==== plotting ====
    # this one will plot all the moments into a single plot
    # any other arguments which will go to imshow
    plot_moment_cmap(cim, ims, line, mom=moment, vmin=0, xlim=[-5,5], ylim=[-8,8], gridmode='pixel')

    # we can output this into a pdf file
    pdfname = 'results/each_line_moment%d.pdf'%moment
    plot_one_at_a_time(cim, ims, line, pdfname, mom=moment)

    # zoom in if you'd like
    pdfname = 'results/each_line_moment%d_z.pdf'%moment
    plot_one_at_a_time(cim, ims, line, pdfname, mom=moment, 
        xlim=[-7, 7], ylim=[-7,7])

    # just here for convenience
    pdfname = 'results/each_line_moment%d_zz.pdf'%moment
    plot_one_at_a_time(cim, ims, line, pdfname, mom=moment,
        xlim=[-3,3], ylim=[-3,3])

if __name__ == '__main__':
    main()

