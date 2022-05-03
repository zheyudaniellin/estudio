""" plot_im2d_masking.py
Plot multiple images while adopting masks for each. 
This is particularly useful for plotting moment 1 and moment 9 images. This sript uses the plot_moment_cmap function from plot_im2d_multi.py
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
import pdb
import etools as eto
from plot_im2d_multi import plot_moment_cmap
au = eto.natconst.au
rad = np.pi / 180.

def plot_one_at_a_time(cim, ims, line, pdfname, **kwargs):
    """ convenient function to iterate through the lines to
    produce one pdf file of continuum with lines one at a time
    """
    pdf = PdfPages(pdfname)

    for i in range(len(ims)):
        fig = plt.figure()
        plot_moment_cmap(cim, [ims[i]], [line[i]], 
            noshow=True, 
            fig=fig, axes=[fig.gca()], 
            gridmode='pixel', **kwargs)
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
#        ['H2CO_3_03-2_02_218.22GHz', 0.5, 'SB'],
        ]

    # system velocity and velocity range for better visualization
    vsys = 5.9
    vwidth = 7.

    # moment for the colormap
    mom = 9
    # 'image' or 'pbcor' for the colormap
    mom_imcor = 'image'

    # moment for the mask
    mask_mom = 8
    # 'image' or 'pbcor' for the mask. It should be 'image'
    mask_imcor = 'image'

    # ==== orgainze the lines ====
    lines = [isrc[0] for isrc in src]
    robust = [isrc[1] for isrc in src]
    baseline = [isrc[2] for isrc in src]

    # ==== reading ====
    # read lines
    ims = []
    for i in range(len(src)):
        # read the color map
        quant = 'mom%d'%mom
        im = eto.image.imageManager()
        fname = eto.kits.get_moment_fits_name(lines[i], baseline[i], robust[i], mom_imcor, mom)
        fname = os.path.join('moments', fname)
        im.read_fits(fname, eto.kits.dpc, quant=quant)
#        setattr(im, '%s_rms'%quant, get_stats(im, quant, [6, 9])['rms'] )

        # read the mask
        quant = 'mom%d'%mask_mom
        fname = eto.kits.get_moment_fits_name(lines[i], baseline[i], robust[i], mask_imcor, mask_mom)
        fname = os.path.join('moments', fname)
        im.read_fits(fname, eto.kits.dpc, quant=quant)
        stats = eto.image.get_stats(im, quant, [6, 9])
        setattr(im, '%s_rms'%quant, stats['rms'] )
        setattr(im, '%s_avg'%quant, stats['avg'] )

        ims.append(im)

    # ==== read continuum ====
    cim = eto.kits.easy_read_continuum('SBLB', 0.5)

    # ==== plotting ====
    pdfname = 'results/each_line_moment%d.pdf'%mom
    plot_one_at_a_time(cim, ims, lines, pdfname,
        mom=mom, mask_mom=mask_mom,
        mthres=1, mground=0.,
        vmin=vsys-vwidth, vmax=vsys+vwidth, cmap=plt.cm.RdBu_r,
        xlim=None,ylim=None)

if __name__ == '__main__':
    main()

