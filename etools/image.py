""" tools to analyze images
- cuts
    1d profiles of an image

- image
    2d image data

adopted units:
- image lengths are in arcseconds
- intensity is in jy/beam
- frequency in Hz
- velocity in km/s
- wavelength in micron
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import pdb
import copy
from astropy.io import fits
from scipy import signal
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline
from . import natconst

class imageCut(object):
    """ an object to record the info of 
    a spatial cut through an image
    """
    def __init__(self):
        self.quantity = []

    def set_quantity(self, quantname, quant):
        """ sets a  2d quantity. the quantity can be anything
        as long as it's in nl by nw
        Parameters
        ----------
        quant : 2d np.ndarray
            the array
        quantname : str
            name of the quantity
        """
        self.quantity.append(quantname)
        setattr(self, quantname, quant)

    def set_laxis(self, laxis):
        self.laxis = laxis

    def set_xy(self, xpnt, ypnt):
        """ xpnt, ypnt in cm
        """
        self.xpnt = xpnt
        self.ypnt = ypnt

    def get_wavelength(self):
        self.w = natconst.cc * 1e4 / self.f

    def set_frequency(self, f):
        self.f = f
        self.nf = len(f)

    def set_dpc(self, dpc):
        self.dpc = dpc

    def set_restfreq(self, restfreq):
        """ set the rest frequency for a line
        """
        self.restfreq = restfreq

    def get_velocity(self):
        """ get the line of sight velocity for lines
        in cm/s
        """
        if hasattr(self, 'f') is False:
            self.get_frequency()

        if hasattr(self, 'restfreq') is False:
            raise ValueError('restfreq must be known to get the velocity')

        self.v = natconst.cc * (self.restfreq - self.f) / self.restfreq / 1e5

    def get_tb(self):
        """ calculate brightness temperature. always in Kelvin 
        """
        self.tb = np.zeros_like(self.I)
        for ii in range(len(self.w)):
            if self.stokes_unit == 'cgs':
                self.tb[...,ii] = cgs_to_tb(self.I[...,ii], self.w[ii])

            elif self.stokes_unit == 'jyppix':
                self.tb[...,ii] = jyppix_to_tb(
                    self.I[...,ii], self.pixel_area, self.w[ii]
                )

            elif self.stokes_unit == 'jypbeam':
                self.tb[...,ii] = jypbeam_to_tb(
                    self.I[...,ii],
                    self.psfarg['bmaj'], self.psfarg['bmin'],
                    self.w[ii]
                )
            else:
                raise ValueError('Internal error. stokes_unit unknown: %s'%self.stokes_unit)

        self.quantity.append('tb')

# =======================================
# image 
# =======================================
class gridManager(object):
    def __init__(self):
        pass

    def setup_from_fits(self, hdr, xdir='west'):
        """ setup the image grid based on the hdr
        xdir : str
            determine if x is to the 'east' or 'west'
        """
        self.nx = hdr['naxis1']
        self.ra = (np.arange(self.nx) - (hdr['crpix1']-1)) * hdr['cdelt1'] + hdr['crval1']

        if xdir == 'east':
            raise ValueError('do not use xdir=east')
            self.dx = hdr['cdelt1'] * 3600.
            self.x = (np.arange(self.nx) - (hdr['crpix1']-1)) * self.dx
        elif xdir == 'west':
            self.dx = - hdr['cdelt1'] * 3600.
            self.x = (np.arange(self.nx) - (hdr['crpix1']-1)) * self.dx
        else:
            raise ValueError('xdir unknown: %s'%xdir)

        self.ny = hdr['naxis2']
        self.dy = hdr['cdelt2'] * 3600.
        self.y = (np.arange(self.ny) - (hdr['crpix2']-1)) * self.dy
        self.dec = (np.arange(self.ny) - (hdr['crpix2']-1)) * hdr['cdelt2'] +hdr['crval2']

        # see if there is a frequency grid
        try:
            self.nf = hdr['naxis3']
            self.f = (np.arange(self.nf) - (hdr['crpix3'] - 1) ) * hdr['cdelt3']+ hdr['crval3']
        except:
            self.nf = 1

        self.fits_hdr = hdr

    def set_wavelength(self, w):
        self.w = w
        self.nf = len(self.w)

    def set_frequency(self, f):
        self.f = f
        self.nf = len(f)

    def get_wavelength(self):
        self.w = natconst.cc * 1e4 / self.f

    def get_frequency(self):
        self.f = natconst.cc * 1e4 / self.w
        self.nf = len(self.f)

    def set_dpc(self, dpc):
        """ distance in parsec
        """
        self.dpc = dpc

    def set_restfreq(self, restfreq):
        """ set the rest frequency for a line
        """
        self.restfreq = restfreq

    def get_pixel_area(self):
        self.pixel_area = self.dx * self.dy

    def get_velocity(self):
        """ get the line of sight velocity for lines
        in km/s
        """
        if hasattr(self, 'f') is False:
            self.get_freqeuncy()

        if hasattr(self, 'restfreq') is False:
            raise ValueError('restfreq must be known to get the velocity')

        self.v = natconst.cc * (self.restfreq - self.f) / self.restfreq / 1e5

    def recenter(self, x0, y0):
        """ reset the center
        x0, y0 : float 
            the center in the current coordinate system 
        """
        self.x -= x0
        self.y -= y0

    def get_mesh(self):
        self.xx, self.yy = np.meshgrid(self.x, self.y, indexing='ij')
        self.rr = np.sqrt(self.xx**2 + self.yy**2)
        self.pp = np.arctan2(self.yy, self.xx)

    def rotate_in_sky(self, phi):
        """ rotate the grid in the sky 
        """
        x1 = self.xx * np.cos(phi) - self.yy * np.sin(phi)
        y1 = self.xx * np.sin(phi) + self.yy * np.cos(phi)
        self.xx, self.yy = x1, y1
        self.pp = np.arctan2(self.yy, self.xx)

class imageManager(object):
    """ 
    Attributes
    ----------
    quantity : list of str
        The list of all the quantities that are set for this imageManager. All the quantities should be in x by y by frequency. 
        Common quantities include: 'I' for continuum or channel map, 'mom0' for moment 0, 'tb' for brightness temperature

    """
    def __init__(self):
        self.quantity = [] # keep track of the physical quantity names

    def set_quantity(self, quantname, quant):
        """ sets a  3d quantity. the quantity can be anything
        as long as it's in nx by ny by nw
        Parameters
        ----------
        quant : 3d np.ndarray
            the array
        quantname : str
            name of the quantity
        """
        self.quantity.append(quantname)
        setattr(self, quantname, quant)

    def trim(self, par):
        """ take out a section of the image
        """
        # trim the x-axis
        xreg = (min(par['xlim']) <= self.grid.x) & (self.grid.x <= max(par['xlim']))
        self.grid.x = self.grid.x[xreg]
        self.grid.nx = len(self.grid.x)
        self.grid.ra = self.grid.ra[xreg]

        # trim the y-axis
        yreg = (min(par['ylim']) <= self.grid.y) & (self.grid.y <= max(par['ylim']))
        self.grid.y = self.grid.y[yreg]
        self.grid.ny = len(self.grid.y)
        self.grid.dec = self.grid.dec[yreg]

        # trim the image
        for q in self.quantity:
            dum = getattr(self, q)[xreg,:,:]
            setattr(self, q, dum[:,yreg,:])

    def get_peak_loc(self, quant='I', iwav=0):
        """ get the location of the peak 
        """
        if quant not in self.quantity:
            raise ValueError('quant unknown: %s'%quant)

        im = getattr(self, quant)[:,:,iwav]

        inx = np.unravel_index(np.nanargmax(im), im.shape)

        return self.grid.x[inx[0]], self.grid.y[inx[1]]

    def get_fn_quant_f(self, quant='I', iwav=0, **kwargs):
        """ returns the interp2d object for just one frequency
        Given the number of pixel required for one image, it's typically impossible to use interp2d and use the meshgrid of x and y. 
        We have to stick to interpolation on a regular grid
        """
        if quant not in self.quantity:
            raise ValueError('the property to be interpolated does not exist: %s'%quant)

        # note that the xaxis is in decreasing ra, which RectBivariateSpline won't allow
        # so we need to shift it
        im2d = getattr(self, quant)[:,:,iwav]
        x = self.grid.x
        y = self.grid.y
        if self.grid.dx < 0:
            x = x[::-1]
            im2d = im2d[::-1,:]
        if self.grid.dy < 0:
            y = y[::-1]
            im2d = im2d[:,::-1]

        fn = RectBivariateSpline(x, y, im2d, **kwargs)

        return fn

    def get_fn_quant(self, quant='I', **kwargs):
        """ prepare interp2d object for each wavelength and keep as attribute
        """
        if quant not in self.quantity:
            raise ValueError('the property to be interpolated does not exist: %s'%quant)

        ff = []
        for ii in range(self.grid.nf):
            fn = self.get_fn_quant_f(quant=quant, iwav=ii, **kwargs)
            ff.append(fn)

        setattr(self, 'fn_%s'%quant, ff)

    def cut(self, laxis, quant=['I'], track='linear',
        trackkw={'x0':0, 'y0':0, 'theta':0},
        mode='fn',
        ):
        """ 
        interpolate a cut across some location 
        Parameters
        ----------
        quant : list of str
            quantity to be sliced
        laxis : 1d ndarray
            some parameter

        track : str
            linear = a linear cut
            ellipse = an elliptical cut

        trackkw : dict
            keyword arguments for the track of the cut 

        mode : str
            'fn' = use the fn_quant method to iterate each wavelength
            'lin' = search for nearest indices and linear interpolate. easier for interpolating throughout wavelength at once

        Returns
        -------
        imageCut object 
        """
        # ==== calculate path ====
        # always in arcsec
        if track == 'input':
            # use arbitrary input
            xpnt = trackkw['x']
            ypnt = trackkw['y']
            # must be the same length as laxis
            if (len(laxis) != len(xpnt)) | (len(laxis) != len(ypnt)):
                raise ValueError("length of x and y must equal laxis when track = 'input'")
        elif track == 'linear':
            """ a linear track
            theta : float, angle east-of-north
            """
            # laxis is the length
            xpnt = trackkw['x0'] + laxis * np.cos(trackkw['theta'] + np.pi/2.)
            ypnt = trackkw['y0'] + laxis * np.sin(trackkw['theta'] + np.pi/2.)

        elif track == 'ellipse':
            # laxis is the phi
            xpnt = trackkw['r'] * np.cos(laxis)
            ypnt = trackkw['r'] * np.sin(laxis)

            # rotate by inclination 
            ypnt, zpnt = rot_inc(trackkw['inc'], ypnt, 0)

            # rotate by some angle
            xpnt, ypnt = rot_eon(trackkw['ang'], xpnt, ypnt)

            xpnt += trackkw['x0']
            ypnt += trackkw['y0']

        else:
            raise ValueError('track unknown')

        # prepare the imageCut object
        cut = imageCut()
        cut.set_laxis(laxis)
        cut.set_xy(xpnt, ypnt)
        cut.set_frequency(self.grid.f)
        cut.get_wavelength()
        cut.stokes_unit = self.stokes_unit

        # ==== interpolation ====
        if mode == 'fn':
            for iquant in quant:
                # check if there is already an interpolation object
                # setup the interpolation if it doesn't have one
                if not hasattr(self, 'fn_%s'%iquant):
                    self.get_fn_quant(quant=iquant)

                fn = getattr(self, 'fn_%s'%iquant)

                # interpolate
                prof = np.zeros([len(laxis), self.grid.nf])
                for i in range(self.grid.nf):
                    prof[:,i] = fn[i](xpnt, ypnt, grid=False)

                cut.set_quantity(iquant, prof)

                del fn, prof

        elif mode == 'lin':
            # find nearest index and interpolate
            # it seems like the fn method is pretty fast already
            pass
        else:
            raise ValueError('cut mode unknown: %s'%mode)

        # some additional attributes
        if hasattr(self.grid, 'dpc'):
            cut.set_dpc(self.grid.dpc)

        if hasattr(self.grid, 'restfreq'):
            cut.set_restfreq(self.grid.restfreq)
            cut.get_velocity()

        return cut

    def set_stokes_unit(self, stokes_unit):
        """ determine the unit of the stokes value
        """
        if stokes_unit not in ['cgs', 'jyppix', 'jypbeam']:
            raise ValueError('stokes_unit unknown: %s'%stokes_unit)
        self.stokes_unit = stokes_unit

    def set_rms(self, rms):
        """ set the rms value
        rms : float
        """
        self.rms = rms

    def read_fits(self, fname, dpc, quant='I'):
        """ read the fits
        """
        hdul = fits.open(fname)
        hdr = hdul[0].header

        # setup the grid
        self.grid = gridManager()
        self.grid.dpc = dpc
        self.grid.setup_from_fits(hdr)

        dat = hdul[0].data.T

        dim = dat.shape
        if len(dim) == 2:
            dat3d = dat[:,:,None]
        elif len(dim) == 3:
            dat3d = dat
        elif len(dim) == 4:
            dat3d = dat[:,:,:,0]
        else:
            raise ValueError('unexpected dimensions: %d'%len(dim))

        if quant == 'I':
            self.set_quantity(quant, dat3d)
            self.set_stokes_unit('jypbeam')
        else:
            self.set_quantity(quant, dat3d)

        self.psfarg = {
                'bmaj':hdr['bmaj']*3600,
                'bmin':hdr['bmin']*3600,
                'bpa':hdr['bpa'],
                'type': 'gaussian',
            }

        hdul.close()

    def get_tb(self):
        """ calculate brightness temperature. always in Kelvin 
        """
        self.quantity.append('tb')

        self.tb = np.zeros_like(self.I)
        for ii in range(len(self.grid.w)):
            if self.stokes_unit == 'cgs':
                self.tb[...,ii] = cgs_to_tb(self.I[...,ii], self.grid.w[ii])

            elif self.stokes_unit == 'jyppix':
                self.tb[...,ii] = jyppix_to_tb(
                    self.I[...,ii], self.grid.pixel_area, self.grid.w[ii]
                )

            elif self.stokes_unit == 'jypbeam':
                self.tb[...,ii] = jypbeam_to_tb(
                    self.I[...,ii],
                    self.psfarg['bmaj'], self.psfarg['bmin'],
                    self.grid.w[ii]
                )
            else:
                raise ValueError('Internal error. stokes_unit unknown: %s'%self.stokes_unit)

    def convert_stokes_unit(self, newunit):
        """ convert the unit for the stokes parameters
        it doesn't do anything if the current unit and the new unit are the same
        newunit : str
        """
        # check if the new unit is the same as the current one
        if self.stokes_unit == newunit:
            return

        # ==== check for necessary arguments ====
        # jyppix
        if (self.stokes_unit == 'jyppix') | (newunit == 'jyppix'):
            if hasattr(self.grid, 'dpc') is False:
                raise ValueError('dpc unknown')

            if hasattr(self.grid, 'pixel_area') is False:
                self.grid.get_pixel_area()

        # jyppbeam
        if (self.stokes_unit == 'jypbeam') | (newunit == 'jypbeam'):
            if hasattr(self, 'psfarg') is False:
                raise ValueError('psfarg unknown')

        # ==== convert the stokes units ====
        if self.stokes_unit == 'cgs':
            if newunit == 'jyppix':
                fac = self.grid.pixel_area / (self.grid.dpc * natconst.pc)**2 / natconst.jy
                if type(fac) == np.ndarray:
                    fac = fac[:,:,None]

                for iquant in ['I','Q','U', 'V', 'pi']:
                    if iquant in self.quantity:
                        setattr(self, iquant, getattr(self, iquant) * fac)

            elif newunit == 'jypbeam':
                # since beam is wavelength dependent
                fac = (self.psfarg['bmaj'] / 3600*natconst.rad) * (self.psfarg['bmin'] / 3600* natconst.rad) * np.pi/4./np.log(2.) / natconst.jy
                for iquant in ['I','Q','U', 'V', 'pi']:
                    if iquant in self.quantity:
                        setattr(self, iquant, getattr(self, iquant) * fac[None,None,:])
            else:
                raise ValueError('newunit unknown: %s'%newunit)

        elif self.stokes_unit == 'jyppix':
            if newunit == 'cgs':
                fac = self.grid.pixel_area / (self.grid.dpc * natconst.pc)**2 / natconst.jy
                if type(fac) == np.ndarray:
                    fac = fac[:,:,None]

                for iquant in ['I']:
                    if iquant in self.quantity:
                        setattr(self, iquant, getattr(self, iquant) / fac)
            elif newunit == 'jypbeam':
                raise ValueError('not done yet')
            else:
                raise ValueError('newunit unknown: %s'%newunit)

        elif self.stokes_unit == 'jypbeam':
            if newunit == 'cgs':
                # since beam is wavelength dependent
                fac = (self.psfarg['bmaj'] / 3600*natconst.rad) * (self.psfarg['bmin'] / 3600*natconst.rad) * np.pi/4./np.log(2.) / natconst.jy
                for iquant in ['I']:
                    if iquant in self.quantity:
                        setattr(self, iquant, getattr(self, iquant) / fac[None,None,:] )
            elif newunit ==  'jyppix':
                raise ValueError('not done yet')
            else:
                raise ValueError('newunit unknown: %s'%newunit)
        else:
            raise ValueError('inner inconsistency. stokes_unit unknown: %s'%self.stokes_unit)

        self.set_stokes_unit(newunit)

    def get_total_flux(self):
        """ total flux as a function of wavelength 
        """
        if hasattr(self.grid, 'pixel_area') is False:
            self.grid.get_pixel_area()
        self.total_flux = np.sum(self.I * self.grid.pixel_area, axis=(0,1))

    def average_channel(self, navg):
        """ average the images across frequency
        navg : int
            the desired number of channels to average
        """
        if navg <= 1:
            return 

        nfull = self.grid.nf // navg
        nres = self.grid.nf -  nfull * navg

        if nres == 0:
            nf = nfull
        else:
            nf = nfull + 1

        newim = np.zeros([self.grid.nx, self.grid.ny, nf])
        newf = np.zeros([nf])

        # Begin iteration 
        for i in range(nfull):
            # find the indices
            inx = np.arange(navg, dtype=int) + i * navg

            newim[:,:,i] = np.mean(self.I[:,:,inx], axis=2)
            newf[i] = np.mean(self.grid.f[inx])

        # Take care of the residual images
        if nres != 0:
            newim[:,:,-1] = np.mean(self.I[:,:,nfull*navg:], axis=2)
            newf[-1] = np.mean(self.grid.f[nfull*navg:])

        # record the changes
        self.I = newim
        self.grid.set_frequency(newf)
        self.grid.get_wavelength()
        self.grid.get_velocity()

    def shift_ground(self, quant, ground):
        """ subtract the quantity by some ground
        """
        if quant not in self.quantity:
            raise ValueError('quantity not found: %s'%quant)

        setattr(self, quant, getattr(self, quant) - ground)

    # ==== plotting ====
    def plot(self, quant='I', iwav=0, ax=None, **kwargs):
        """ plot the 2d quantity
        """
        if ax is None:
            ax = plt.gca()
        im2d = getattr(self, quant)[:,:,iwav]
        pc = ax.pcolormesh(self.grid.x, self.grid.y, im2d.T, **kwargs)

    def plot_beam(self, ax=None, beamxy=None, axis_unit='arcsec', rotsky=0, facecolor='w'):
        """ plot the gaussian beam size
        Parameters
        ----------
        implot : ax map? 
        beamxy : tuple
            location of the center of the beam 
        axis_unit : string
            unit the beam should be and also for beamxy
        rotsky : float
            angle of rotation in the sky in radians. angle is counter-clockwise from the x-axis
        """
        if self.psfarg['type'] != 'gaussian':
            raise ValueError('psf type not available')

        if ax is None:
            ax = plt.gca()

        # determine size of ellipse to overplot
        # conversion factor from arcsec 
        if axis_unit == 'cm':
            fac = self.grid.dpc * natconst.au
        elif axis_unit == 'au':
            fac = self.grid.dpc
        elif axis_unit == 'arcsec':
            fac = 1.
        else:
            raise ValueError('axis_unit unknown: %s'%axis_unit)

        ewidth = self.psfarg['bmin'] * fac
        eheight = self.psfarg['bmaj'] * fac

        # center of the beam
        if beamxy is None:
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            ex = xmin + 0.8 * (xmax - xmin)
            ey = ymin + 0.2 * (ymax - ymin)
        else:
            ex = beamxy[0]
            ey = beamxy[1]

        # angle
#        ang = 180. - self.psfarg['bpa'] - rotsky / natconst.rad
        ang = self.psfarg['bpa']

        # set up the ellipse
        ells = Ellipse( xy=(ex,ey),
            width=ewidth, height=eheight, angle=ang)
        ells.set_facecolor(facecolor)
        ells.set_fill(True)
        ax.add_patch(ells)

# ======
# other plotting tools
# ======
def plot_length(ax, x, y, length, txt, color='w', **kwargs):
    """ plot a scale length in a plot
    Parameters
    ----------
    ax : matplotlib axes object
    x, y : float
        The center of the line in units corresponding to the current plot
    length : float
        The length of the line in units corresponding to the current plot
    txt : str
        The text used to annotate the line
    """
    pltx = x + np.array([-1,1]) * 0.5 * length
    plty = y + np.array([0,0])
    ax.plot(pltx, plty, color=color, **kwargs)
    ax.text(x, y, txt, va='bottom', ha='center', color=color)

def plot_beam(ax, x, y, bmaj, bmin, bpa, facecolor='w'):
    """ plot the beam
    ax : matplotlib axes object
    x, y : float
        the location of the ellipse
    bmaj, bmin : float
        The beam major and minor axis length in units corresponding to the current plot
    bpa : float
        The beam angle east-of-north
    facecolor : str
        The color of the beam 
    """
    ang = bpa
    ells = Ellipse(xy=(x,y), width=bmin, height=bmaj, angle=ang)
    ells.set_facecolor(facecolor)
    ells.set_fill(True)
    ax.add_patch(ells)

# ====================
# image unit conversions
# ====================
def calc_beam_solid_angle(bmaj, bmin):
    """ bmaj, bmin in arcsec
    """
    return (bmaj / 3600. * natconst.rad) * (bmin / 3600. * natconst.rad) * np.pi/4. / np.log(2.)

def jyppix_to_cgs(jyppix, dxy):
    """
    convert jy/pixel image to inu in cgs
    Parameters
    ----------
    jyppix       : ndarray
        the image in jy/pixel units
    dxy     : float
        the solid angle size the pixel in arcsec**2
    """
    solid_angle = dxy / 3600.**2 * natconst.rad**2
    return jyppix * natconst.jy / solid_angle

def cgs_to_jyppix(cgsim, dxy):
    solid_angle = dxy / 3600.**2 * natconst.rad**2
    return cgsim * solid_angle / natconst.jy

def jypbeam_to_cgs(jypbeam, bmaj, bmin):
    """ convert image of jy/beam to cgs units
    """
    beam = calc_beam_solid_angle(bmaj, bmin)
    return jypbeam * natconst.jy / beam

def cgs_to_jypbeam(cgsim, bmaj, bmin):
    beam = calc_beam_solid_angle(bmaj, bmin)
    return cgsim * beam / natconst.jy

def cgs_to_tb(cgsim, wav_micron):
    """ convert image in cgs units to brightness temperature
    """
    freq = natconst.cc * 1e4 / wav_micron
    ld2 = (wav_micron*1e-4)**2
    hnu = natconst.hh * freq
    hnu3_c2 = natconst.hh * freq**3 / natconst.cc**2
    tb = hnu / natconst.kk / np.log(2. * hnu3_c2 / abs(cgsim) + 1.) * np.sign(cgsim)
    return tb

def tb_to_cgs(tb, wav_micron):
    """ convert image in brightness temperature to cgs units
    basically the planck function
    """
    freq = natconst.cc * 1e4 / wav_micron
    ld2 = (wav_micron * 1e-4)**2
    hnu = natconst.hh * freq
    hnu3_c2 = natconst.hh * freq**3 / natconst.cc**2
    cgs = 2 * hnu3_c2 / (np.exp(hnu / natconst.kk / abs(tb)) - 1.) * np.sign(tb)
    return cgs

def jyppix_to_tb(jyppix, dxy, wav_micron):
    """
    converts jy/pixel image to brightness temperature 
    """
    cgsim = jyppix_to_cgs(jyppix, dxy)
    tb = cgs_to_tb(cgsim, wav_micron)
    return tb

def jypbeam_to_tb(jypbeam, bmaj, bmin, wav_micron):
    """
    converts jy/beam image to brightness temperature 
    """
    cgsim = jypbeam_to_cgs(jypbeam, bmaj, bmin)
    tb = cgs_to_tb(cgsim, wav_micron)
    return tb

def tb_to_jypbeam(tb, bmaj, bmin, wav_micron):
    """ convert brightness temperature to jy/beam
    """
    cgsim = tb_to_cgs(tb, wav_micron)
    jypbeam = cgs_to_jypbeam(cgsim, bmaj, bmin)
    return jypbeam

# ====================
# basic image calculations
# ====================
def get_stats(im, quant, rlim, iwav=None):
    """ calculate the rms based on some region
    im : imageManager object
    quant : str
        name of the desired quantity, like 'I', 'mom0'
    rlim : list of 2 floats
        the range in radius in arcsec
    iwav : int
        the index in wavelength. If None, then the stats will be across wavelength
    """
    xx, yy = np.meshgrid(im.grid.x, im.grid.y, indexing='ij')
    rr = np.sqrt(xx**2 + yy**2)
    reg = (rlim[0] <= rr) & (rr <= rlim[1])
    if iwav is None:
        vals = getattr(im, quant)[reg,:]
    else:
        vals = getattr(im, quant)[reg,iwav]

    return {'rms':np.nanstd(vals), 'avg':np.nanmean(vals)}

