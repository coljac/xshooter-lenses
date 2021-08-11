#!env python
import re
import os
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.io.fits import Header

def get_fitsfiles(directory, searchstring):
    files = []
    for filter in "irg":
        ffind = re.compile(r"_"+filter+"[._]")
        for filen in os.listdir(directory):
            if searchstring in filen and ".fits" in filen and ffind.search(filen):
                files.append(filen)

    iband, rband, gband = files
    return files

# Fitsfile OR header object
def pix_to_sky(fitsfile, x, y):
    w = WCS(fitsfile)
    ra, dec = w.all_pix2world(x, y, 0)
    return ra, dec

def sky_to_pix(fitsfile, ra, dec):
    w = WCS(fitsfile)
    x, y = w.all_world2pix(ra, dec, 1)
    return int(x), int(y)

def get_header(fitsfile):
    ff = fits.open(fitsfile)
    return ff[0].header

def header_from_string(headerstring):
    return Header.fromstring(headerstring)

# def pointing_to_fits(pointing, filterband="g", directory="."):
    # for filen in os.listdir(directory):
        # if fnmatch.fnmatch(filen, '*' + "_" + filterband + "_" + pointing + '*.fits'):
            # return filen
    # return None

def lens_to_coords(ra, dec, pointing):
    fitsname = pointing_to_fits(pointing)
    if not fitsname:
        return None
    return sky_to_pix(fitsname, ra, dec)

def pointing_pixels_to_coords(x, y, pointing):
    fitsname = pointing_to_fits(pointing)
    if not fitsname:
        return None
    return pix_to_sky(fitsname, x, y)

def make_cutout(infile, x, y, size, outfile):
    startx, starty, width, height = x-size/2, y-size/2, size, size

    hdulist = fits.open(infile, memmap=True)
    headers = hdulist[0].header
    indata = hdulist[0].data

    outdata = indata[starty:(starty+height), startx:(startx+width)]

    head = headers.copy()
    outfits = fits.PrimaryHDU(data=outdata, header=head)
    outfits.writeto(outfile, clobber=True)

    hdulist.close()

