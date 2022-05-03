# script_moments.py
# casa script to calculate the moments
import os

### the naming scheme
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

def get_imagename(prefix, src, imcor, baseline='SB'):
    """ get the image name
    baseline : str
        the combo of baseline
        'SB', 'SBLB', 'LB'
    """
    namelist = []
    for ikey in src.keys():
        for irobust in src[ikey]:
            imagename = '%s'%prefix + '_%s'%baseline + '_%s'%ikey + '_robust_%.1f'%irobust + '.%s'%imcor
            namelist.append(imagename)

    return namelist

def iterate_moments(moments, imagenames, datdir, outdir):
    """ iterate through each image to calculate the moments
    Note that if moments is a list of only one value, then casa will not give the appropriate suffix, so this will do it for you. 
    This will only keep the fits file and the casa image will be deleted
    """
    moment_tags = [moment_base[str(i)] for i in moments]

    for iname in imagenames:
        # determine the image name
        image = os.path.join(datdir, iname + '.fits')

        # determine output prefix
        outprefix = os.path.join(outdir, iname)

        # remove any pre-existing files
        for itag in moment_tags:
            fname = outprefix + '.' + itag
            os.system('rm -rf ' + fname)
            fname = outprefix + '.' + itag + '.fits'
            os.system('rm -rf ' + fname)

        # if give the appropriate suffix if only one element
        # since casa doesn't do it
        if len(moments) == 1:
            outfile = outprefix + '.' + moment_tags[0]
        else:
            outfile = outprefix

        # call moments
        immoments(imagename=image, moments=moments, outfile=outfile)

        # export fits
        for i, imoment in enumerate(moments):
            image = outprefix + '.' + moment_tags[i]
            fitsimage = image + '.fits'
            exportfits(imagename=image, fitsimage=fitsimage)

            # delete the casa image. we only want the fits file
            os.system('rm -rf ' + image)

# ==== settings ====
baseline = 'SB'
prefix = 'IRAS04302'

# this is the directory where the exported data is
datdir = '/lustre/cv/projects/edisk/IRAS04302/export_%s'%baseline

# the directory for the output 
outdir = 'moments'

# change these based on what you have
# each entry is the name of the line for the fits and the robust parameters
# comment or uncomment based on what you need 
if baseline == 'SB':
    src = {
        '12CO':[0.5],
        '13CO':[-0.5, 0.5],
        'C18O':[-0.5, 0, 0.5], 
        'SO':[0.5], 
#        'H2CO_3_21-2_20_218.76GHz':[0.5], 
#        'H2CO_3_03-2_02_218.22GHz':[0.5], 
#        'H2CO_3_22-2_21_218.47GHz':[0.5], 
#        'c-C3H2_217.82':[0.5], 
#        'cC3H2_217.94':[0.5],
#        'cC3H2_218.16':[0.5], 
#        'DCN':[0.5], 
#        'CH3OH':[0.5],
#        'SiO':[0.5], 
        }
else:
    src = {
        '12CO':[0.5],
        '13CO':[0.5],
        'C18O':[0.5], 
        'SO':[0.5], 
        'H2CO_3_21-2_20_218.76GHz':[0.5], 
        'H2CO_3_03-2_02_218.22GHz':[0.5], 
        'H2CO_3_22-2_21_218.47GHz':[0.5], 
        'c-C3H2_217.82':[0.5], 
        'cC3H2_217.94':[0.5],
        'cC3H2_218.16':[0.5], 
        'DCN':[0.5], 
        'CH3OH':[0.5],
        'SiO':[0.5], 
        }

# choose either 'image' or 'pbcor'
imcor = 'image'

# this can be a list of any of the moments you want, e.g. [0,1,2,8,9,10]
moments = [9]

# ==== calculation ====
namelist = get_imagename(prefix, src, imcor, baseline=baseline)
iterate_moments(moments, namelist, datdir, outdir)


