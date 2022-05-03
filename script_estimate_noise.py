""" casa script to estimate noise for a series of images
Estimating the noise in advance can speed up future plotting.

The results will be written in the 'results' directory. 
This is basically the same as in the reduction_utils3.py from edisk

To use this script, enter casa and type
execfile('script_estimate_noise.py')

Be sure to check out the functions get_continuum() and get_lines()
as there are certain parameters to change based on each source
"""
import numpy as np
import pickle
import pdb
import os

def estimate_stats(imagename, disk_mask, noise_mask, noisechan=''):
    """ based on estimate_SNR from reduction_utils3.py
    calculates the following from an image or image cube
    - beam major and minor axis, beam pa
    - total flux within the disk 
    - peak intensity within the disk
    - rms in the noise regions
    - signal to noise ratio of the peak intensity
    The function prints the results and also stores it to a dictionary

    Parameters
    ----------
    imagename : str
    disk_mask : str
        in CASA region format
    noise_mask : str

    Returns
    -------
    dictionary
    """
    # beam info
    headerlist = imhead(imagename, mode = 'list')
    beammajor = headerlist['beammajor']['value']
    beamminor = headerlist['beamminor']['value']
    beampa = headerlist['beampa']['value']
    print("#====")
    print("#%s" % imagename)
    print("#Beam %.3f arcsec x %.3f arcsec (%.2f deg)" % (beammajor, beamminor, beampa))

    # properties within disk mask
    disk_stats = imstat(imagename = imagename, region = disk_mask)
    disk_flux = disk_stats['flux'][0]
    print("#Flux inside disk mask: %.2f mJy" % (disk_flux*1000,))
    peak_intensity = disk_stats['max'][0]
    print("#Peak intensity of source: %.2f mJy/beam" % (peak_intensity*1000,))

    # noise related
    noise_stats = imstat(imagename=imagename, region=noise_mask, chans=noisechan)
    rms = noise_stats['rms'][0]
    print("#rms: %.2e mJy/beam" % (rms*1000,))

    # signal to noise of peak 
    SNR = peak_intensity / rms

    print('#snr: %.2e' % SNR)

    return {'bmaj':beammajor, 'bmin':beamminor, 'bpa':beampa, 
        'flux':disk_flux, 'peak':peak_intensity, 
        'rms':rms, 'snr':SNR, }

# ==== continuum ====
def get_continuum(mask_ra, mask_dec, prefix, baseline, datdir):
    """ read the continuum and estimate the noise
    baseline : str
    """
    mask_pa  = 175.5        # position angle of mask in degrees
    mask_maj = 3.0  # semimajor axis of mask in arcsec
    mask_min = 1.0  # semiminor axis of mask in arcsec

    common_mask = 'ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]' % \
              (mask_ra, mask_dec, mask_maj, mask_min, mask_pa)
    """ Define a noise annulus, measure the peak SNR in map """
    noise_annulus = "annulus[[%s, %s],['%.2farcsec', '8.0arcsec']]" % \
                (mask_ra, mask_dec, 2.0*mask_maj)

    stats = {}
    nameslist = []
    for robust in [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]:
        keyname = prefix + '_%s'%baseline + '_continuum_robust_' + str(robust) + '.image.tt0'
        imagename = os.path.join(datdir, 'export_%s'%baseline, keyname + '.fits')

        istat = estimate_stats(imagename, common_mask, noise_annulus)
        stats[keyname] = istat
        nameslist.append(keyname)

    if baseline == 'SBLB':
        for taper in ['1000klambda', '2000klambda', '3000klambda']:
            for robust in [1.0, 2.0]:
                keyname = prefix+'_%s'%baseline + '_continuum_robust_'+str(robust)+'_taper_'+taper + '.image.tt0'
                imagename = os.path.join(datdir, 'export_%s'%baseline, keyname + '.fits')
 
                istat = estimate_stats(imagename, common_mask, noise_annulus)

                stats[keyname] = istat
                nameslist.append(keyname)

    stats['disk_mask'] = common_mask
    stats['noise_mask'] = noise_annulus
    stats['nameslist'] = nameslist

    return stats

# ==== lines ====
def get_lines(mask_ra, mask_dec, prefix, baseline, datdir):
    """ calculate the basic stats for lines
    """
    mask_pa  = 175.5        # position angle of mask in degrees
    mask_maj = 6.  # semimajor axis of mask in arcsec
    mask_min = 6.0  # semiminor axis of mask in arcsec

    common_mask = 'ellipse[[%s, %s], [%.1farcsec, %.1farcsec], %.1fdeg]' % \
              (mask_ra, mask_dec, mask_maj, mask_min, mask_pa)

    """ Define a noise annulus, measure the peak SNR in ma """
    # usually the first few channels are completely emissionless, 
    # but I'm still using a mask just in case
    rmin = '3.0arcsec'
    rmax = '8.0arcsec'
    noise_annulus = "annulus[[%s, %s],[%s, %s]]" % \
                (mask_ra, mask_dec, rmin, rmax)
    noisechan = '0~5' # the channels for the noise

    if baseline == 'SB':
        src = {
            'C18O':[-0.5, 0, 0.5], 
            '13CO':[-0.5, 0.5], 
            '12CO':[0.5], 
            'SO':[0.5], 
            'H2CO_3_21-2_20_218.76GHz':[0.5], 
            'H2CO_3_03-2_02_218.22GHz':[0.5], 
            'H2CO_3_22-2_21_218.47GHz':[0.5], 
            'c-C3H2_217.82':[0.5], 
            'cC3H2_217.94':[0.5], 
            'cC3H2_218.16':[0.5],
            'DCN':[0.5], 
            'CH3OH':[0.5],
            'SiO':[0.5]
            }
    else:
        src = {
            'C18O':[0.5],
            '13CO':[0.5],
            '12CO':[0.5],
            'SO':[0.5],
            'H2CO_3_21-2_20_218.76GHz':[0.5],
            'H2CO_3_03-2_02_218.22GHz':[0.5],
            'H2CO_3_22-2_21_218.47GHz':[0.5],
            'c-C3H2_217.82':[0.5],
            'cC3H2_217.94':[0.5],
            'cC3H2_218.16':[0.5],
            'DCN':[0.5],
            'CH3OH':[0.5],
            'SiO':[0.5]
            }

    stats = {}
    nameslist = []
    for ikey in src.keys():
        robust = src[ikey]
        for irobust in robust:
            keyname = prefix + '_%s'%baseline + '_%s'%ikey + '_robust_%.1f'%(irobust) + '.image'
            imagename = os.path.join(datdir, 'export_%s'%baseline, keyname + '.fits')

            istat = estimate_stats(imagename, common_mask, noise_annulus)
            stats[keyname] = istat
            nameslist.append(keyname)

    stats['disk_mask'] = common_mask
    stats['noise_mask'] = noise_annulus
    stats['nameslist'] = nameslist

    return stats

# ==== settings ===
### this is the absolute path to the directory that contains the exported stuff like 'export_SB', 'export_SBLB'
datdir = '/lustre/cv/projects/edisk/IRAS04302'

### the prefix of your source
prefix = 'IRAS04302'

### determine if you want to do 'SB' or 'SBLB'
baseline = 'SBLB'

### do the continuum? 
# this is necessary
do_continuum = True

### do the line? 
# this isn't necessary yet, so you can skip this if you want to
do_lines = False

# ==== determine the center of the masks ====
mask_ra, mask_dec = '04:33:16.49977', '+22.53.20.225224'

# ==== get the continuum results ====
### be sure to look at the get_continuum() function for masking parameters
if do_continuum:
    cont = get_continuum(mask_ra, mask_dec, prefix, baseline, datdir)

    # write to file
    fname = 'results/image_stats_continuum_%s.pickle'%baseline
    with open(fname, 'wb') as f:
        pickle.dump(cont, f, protocol=pickle.HIGHEST_PROTOCOL)

# ==== get the line results ====
### be sure to look at the get_lines() function for masking parameters
if do_lines:
    lines = get_lines(mask_ra, mask_dec, prefix, baseline, datdir)

    # write to file
    fname = 'results/image_stats_line_%s.pickle'%baseline
    with open(fname, 'wb') as f:
        pickle.dump(lines, f, protocol=pickle.HIGHEST_PROTOCOL)


