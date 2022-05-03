""" kits.py
Some basic tool kits that will streamline plotting
"""
import os
import numpy as np
import pickle
import pdb
from . import image

# keep a record of the names of the moments
moment_base = {
    '-1':'average',
    '0':'integrated',
    '1':'weighted_coord',
    '2':'weighted_dispersion_coord',
    '3':'median',
#    '4': # this will not work
    '5':'standard_deviation',
    '6':'rms',
    '7':'abs_mean_dev',
    '8':'maximum',
    '9':'maximum_coord',
    '10':'minimum'
    }

# ==== hard code these since it's the same for each source ====
# the prefix name for your source
prefix = 'IRAS04302'

# the adopted distance
dpc = 140 

# the directory where 'export_SB' and 'export_SBLB' are located
rawdir = '/lustre/cv/projects/edisk/IRAS04302/'

# =============================================================

def get_moment_fits_name(line, baseline, robust, imcor, moment, 
    prefix=prefix):
    """ get the fits name of the moment image
    line : str
    baseline : str
    robust : float
    imcor : str
    moment : int
    """
    fname = '%s_'%prefix + baseline + '_%s'%line + '_robust_%.1f'%robust + '.%s'%imcor

    itag = moment_base['%d'%moment]
    fname += '.' + itag

    fname += '.fits'
    return fname

def easy_read_continuum(baseline, robust, imcor='image', 
    prefix = prefix, 
    dpc = dpc, 
    rawdir = rawdir, 
    ):
    """ read the continuum
    Parameters
    ----------
    baseline : str
    robust : float
    imcor : str
        either 'image' or 'pbcor'
    prefix : str
        The prefix of you source
    dpc : float
        Distance to source in pc
    rawdir : str
        The directory where the 'export_SB' and 'export_SBLB' which are where all the continuum and image cube fits files are kept
    """
    # quick checks
    if baseline not in ['SB', 'SBLB']:
        raise ValueError('baseline unknown')

    if imcor not in ['image', 'pbcor']:
        raise ValueError('image or pbcor unknown')

    # read data
    fname = os.path.join(rawdir, 'export_%s'%baseline, 
        '%s_%s_continuum_robust_%.1f.%s.tt0.fits'%(prefix, baseline, robust, imcor) )
    im = image.imageManager()
    im.read_fits(fname, dpc)

    # read noise
    fname = os.path.join('results', 'image_stats_continuum_%s.pickle'%baseline)
    with open(fname, 'rb') as f:
        stats = pickle.load(f)
    ikey = '%s_%s_continuum_robust_%.1f.image.tt0'%(prefix, baseline, robust)
    im.rms = stats[ikey]['rms']

    return im
