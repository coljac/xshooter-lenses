import matplotlib
matplotlib.use("Agg")
import sys
import os
import fitsutils as fu
import time
import pandas as pd
import numpy as np
import coltools as ct
import astropy.io.fits as fits
from PIL import Image, ImageDraw, ImageOps
import math
import astropy.coordinates as coord
import astropy.units as u
import argparse


from prettytable import PrettyTable
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import rcParams
from photutils import centroid_com, centroid_1dg, centroid_2dg
import datetime
import pickle as pkl
import humvi
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib_scalebar import dimension
import stars as starutils

asec = dimension._Dimension("arcsec", '$^{\\prime}$$^{\\prime}$')
flipit = True
basedir = "."
slitwidth = 5
slitheight = [42]
dims = 100
slitcolors = ["yellow" ,"red"]
# astroquery_imported = False

# def import_astroquery(imported=astroquery_imported):
#     print("IMPORT")
#     if imported:
#         return
# astroquery_imported =True

def region_to_observation(lens, region_file, number=1):
    # def __init__(self, lens, best_slit, nod_slits, region_string, star_location=None, times=None):
    region_text = ct.fal(region_file)
    num = 0
    box = [] # None
    line = None
    pa=None
    for l in region_text:
        if l.startswith("box"):
            # if box is None:
            box.append(l.split("(")[1])
        elif l.startswith("line"):
            num += 1
            if num == 1 and number == 2:
                box = []
                continue
            else:
                line = l.split("(")[1]
                break

    star_location = [float(x) for x in line.split(",")[0:2]]

    slit_location = [float(x) for x in box[0].split(",")[0:2]]
    nod_location = [float(x) for x in box[1].split(",")[0:2]]
    pa = float(box[0].split(")")[0].split(",")[-1])

    pa = 180 - pa  # Needed?

    position = Position((0, 0), pa, radec=slit_location)
    # position.from_ra_dec(lens.header, *slit_location)

    slit = Slit(position, None)
    print(slit_location)
    print(position.get_ra_dec(lens.header))


    print("Pos:", position)
    print("slit:", slit_location)
    print("nod:", nod_location)
    position_nod = Position((0, 0), pa, radec=nod_location)
    position_nod.from_ra_dec(lens.header, *nod_location)
    nod_slit = Slit(position_nod, None)

    print("Read:", pa, slit_location, star_location)
    # def __init__(lens, best_slit, nod_slits, region_string, star_location=None, times=None):
    return Observation(lens, slit, [slit, nod_slit], "", star_location=star_location, times=None)


# def to_png_redo(image_data, save_to, extra_layer=None, band="all", observations=None, crop=True, size=512):
def to_png_redo(lens, starloc, slits, band="g", draw_both=True):
    data = {}
    header = None
    for i, band in enumerate("irg"):
        # data["rgb"[i]] = lens.data_1024(band, flip=True)
        d, h = lens.data_1024(band, flip=False)
        data['rgb'[i]] = fits.PrimaryHDU(d)
        data['rgb'[i]].header = h
        # data["rgb"[i]], header = fits.PrimaryHDU(lens.data_1024(band, flip=True))

        # data = lens.data_1024(band, flip=True).data

    # data['r'] = fits.PrimaryHDU(image_data[2][0])
    # data['g'] = fits.PrimaryHDU(image_data[1][0])
    # data['b'] = fits.PrimaryHDU(image_data[0][0])
    # fits_headers = [x[1] for x in image_data]

    # data['r'].header = fits_headers[2]
    # data['g'].header = fits_headers[1]
    # data['b'].header = fits_headers[0]

    print("Got data:", data[band].shape)

    rscale, gscale, bscale, Q, alpha = (1.0, 1.2, 2.8, 1.0, 0.03)
    masklevel, saturation, offset, backsub, vb = None, 'white', 0.0, False, False

    if band != "all":
        for b in "rgb":
            data[b] = data[band]

    # Sigh. :(
    # if flipit:
    #     for b in "rgb":
    #         data[b].data = np.flipud(data[b].data)

    output_img = humvi.compose(data['r'], data['g'], data['b'], scales=(rscale, gscale, bscale), Q=Q, alpha=alpha,
                               masklevel=masklevel, saturation=saturation, offset=offset, backsub=backsub,
                               vb=vb, outfile=None)
    output_img = ImageOps.invert(output_img)
    draw = ImageDraw.Draw(output_img)


    # First, the star
    size = 1024
    print("Got star loc:", starloc)
    header = data['g'].header
    star_start = Position((0, 0), 0)
    star_start.from_ra_dec(header, *starloc)
    print("Star start:", star_start)

    star_c = star_start.center()[1], size - star_start.center()[0]
    print("Got star centre", lens.name, star_c)

    if star_c[0] < 0 or star_c[1] < 0 or star_c[0] > size or star_c[1] > size:
        print(lens.name, "OFF THE CHARTS!!!", star_c)
        return None

    # draw.rectangle((0, 0, 100, 100), width=3, outline="red")
    # draw.rectangle((0, 0, 500, 500), width=3, outline="blue")
    # draw.rectangle((500, 500, 700, 700), width=3, outline="green")
    # draw.rectangle((463, 319, 500, 400), width=3, outline="purple")

    box_size = 20
    draw.rectangle([(star_c[0] - box_size, star_c[1] - box_size),
         (star_c[0] + box_size, star_c[1] + box_size)],
                     width=2,
                     outline="red")
    colors = ["red", "blue", "red", "blue"]

    slits[0].set_from_ra_dec(header)
    slits[1].set_from_ra_dec(header)

    slit_centre = slits[0] # Pass as positions?
    # slit_centre.set_from_ra_dec(header, slits[0].ra, slits[0].dec)
    # slit_centre.set_from_ra_dec(header)

    line_end = (slit_centre.x, slit_centre.y)
    # line_end = resize_point(*pos_end, lens, to=size)
    draw.line([(star_start.center()[1], size - star_start.center()[0]), line_end], fill="red")

    # Draw slit(s)
    for i, slit in enumerate(slits):
        # draw.rectangle([(slit.x, slit.y), (slit.x+box_size, slit.y + box_size)],
        #                width=2,
        #                outline="blue")
        # Flip angle?
        slit.dims = 1024
        print("Center for image:", slit.center())
        slit.set_from_ra_dec(header)
        print("Center for image:", slit.center())
        # Here's where all the trouble is.
        print("But ra dec is ", slit.get_ra_dec(header), slit.x, slit.y)
        output_img = rotate(slit.theta, slitwidth, slitheight[0],
                            output_img.size[0],
                            rect_center=slit.center(),
                                    # resize_point(*slit.center(),
                                    # lens, to=size),
                            image=output_img, fill=None, outline=colors[i])


    # CROP
    crop_middle = star_c
    crop_value = 512
    print(">>>", crop_middle, line_end)
    sy = crop_middle[0] - crop_value
    sx = crop_middle[1] - crop_value

    arcsec_45 = int(45 / .263)

    draw.rectangle([
        (crop_middle[0] - arcsec_45, crop_middle[1] - arcsec_45),
        (crop_middle[0] + arcsec_45, crop_middle[1] + arcsec_45)
    ], outline="blue")

    output_img = output_img.crop(box=(sy, sx, crop_middle[0] + crop_value, crop_middle[1] + crop_value))

    return output_img

def to_png(image_data, save_to, extra_layer=None, band="all", observations=None, crop=True, size=512):
    data = {}
    data['r'] = fits.PrimaryHDU(image_data[2][0])
    data['g'] = fits.PrimaryHDU(image_data[1][0])
    data['b'] = fits.PrimaryHDU(image_data[0][0])
    fits_headers = [x[1] for x in image_data]

    data['r'].header = fits_headers[2]
    data['g'].header = fits_headers[1]
    data['b'].header = fits_headers[0]

    print("Got data:", data[band].shape)

    rscale, gscale,bscale,Q, alpha =  (1.0, 1.2, 2.8, 1.0, 0.03)
    masklevel, saturation, offset, backsub, vb = None, 'white', 0.0, False, False
    if band != "all":
        for b in "rgb":
            data[b] = data[band]

    # Sigh. :(
    if flipit:
        for b in "rgb":
            data[b].data = np.flipud(data[b].data)

    output_img = humvi.compose(data['r'], data['g'], data['b'], scales=(rscale,gscale,bscale), Q=Q, alpha=alpha, 
      masklevel=masklevel, saturation=saturation, offset=offset, backsub=backsub, 
            vb=vb, outfile=None)
    output_img = ImageOps.invert(output_img)
    draw = ImageDraw.Draw(output_img)
    if observations is not None:
        lens = observations[0].lens
        lensname = lens.name
        starloc = observations[0].star_location
        print("Got star loc:", observations[0].star_location)
        # if lensname == "DESJ0003-3348":
        #     starloc = 0.8320599852,-33.7937746
        #     print("Doing the special one.")

        header = data['g'].header
        star_start = Position((0, 0), 0)
        star_start.from_ra_dec(header, *starloc)
        star_c = star_start.center()[1], size-star_start.center()[0]

        print("Got star centre", observations[0].lens.name, star_c)
        if star_c[0] < 0 or star_c[1] < 0 or star_c[0] > size or star_c[1] > size:
            print(observations[0].lens.name, "OFF THE CHARTS!!!", star_c)
            return None

        draw.rectangle([(star_c[0]-10, star_c[1]-10),
            (star_c[0]+10, star_c[1]+10)],
               outline="red")
        colors = ["red", "blue", "red", "blue"]
        if len(observations) == 1:
            print("Got one observation. with this many slits", len(observations[0].nod_slits))
            print("slit centres: ", observations[0].nod_slits[0].position.center())
            print("slit centres: ", observations[0].nod_slits[1].position.center())
            pos_end = 0.5*(np.array(observations[0].nod_slits[0].position.center()) \
                    + np.array(observations[0].nod_slits[1].position.center()))
            slits = [observations[0].nod_slits[0], observations[0].nod_slits[1]]
        else:
            pos_end = 0.5*(np.array(observations[0].nod_slits[0].position.center()) \
                    + np.array(observations[1].nod_slits[0].position.center()))
            slits = [observations[0].nod_slits[0], observations[0].nod_slits[1],observations[1].nod_slits[0],
observations[1].nod_slits[1]]

        for i, slit in enumerate(slits):
            output_img = rotate(slit.position.theta, slitwidth, slitheight[0],
                                output_img.size[0],
                                rect_center=resize_point(*slit.position.center(),
                                                         lens, to=size),
                                image=output_img, fill=None, outline=colors[i])
        line_end = resize_point(*pos_end, lens, to=size)
        draw.line([(star_start.center()[1], size-star_start.center()[0]), line_end], fill="red")
        print("This: XXXX", star_c)

    crop_middle = star_c
    crop_value = 512
    print(">>>", crop_middle, line_end)
    sy = crop_middle[0] - crop_value
    sx = crop_middle[1] - crop_value 

    arcsec_45 = int(45/.263)

    draw.rectangle([
        (crop_middle[0] - arcsec_45, crop_middle[1] - arcsec_45), 
        (crop_middle[0] + arcsec_45, crop_middle[1] + arcsec_45)
        ],outline="blue")

    output_img = output_img.crop(box=(sy, sx,crop_middle[0]+ crop_value, crop_middle[1]+crop_value))

    if False and crop:
        distance_to_middle = np.abs(512 - np.array(star_c))
        # print(star_c)
        # print(distance_to_middle)
        crop_value = distance_to_middle.max() + 100
        # print(crop_value)
        # min_ = np.array(star_c).min()
        # max_ = size - np.array(star_c).max()
        # print(min_, max_, star_c)
        # crop_value = max(200, min(min_, max_)+50)
        # print(crop_value)
        # crop_value = min(min_, max_) - 20
        # crop_value = min(70, crop_value) 
        # if np.array(star_c).min() > 50 and np.array(star_c).max() > size-50:
            # from_ = 50
            # to = size-50
        # center on lens
        # output_img = output_img.crop(box=(crop_value, crop_value, size-crop_value, size-crop_value))
        # center on star

        crop_middle = star_c
        sy = crop_middle[0] - crop_value
        sx = crop_middle[1] - crop_value 
        output_img = output_img.crop(box=(sy, sx,crop_middle[0]+ crop_value, crop_middle[1]+crop_value))
    return output_img

## This won't work, because the bigger images don't have the smaller ones at the 
def resize_point(y, x, lens, from_=100, to=512):
    # ra, dec = fu.pix_to_sky(lens.header, x, 100-y)

    # ra, dec = fu.pix_to_sky(lens.header, x, y)
    ra, dec = fu.pix_to_sky(lens.header, y, x)

    # ra, dec = fu.pix_to_sky(lens.header, y, x)
    print(x, y, ra, dec, "resized")
    if False and to == 2048:
        header = lens.data_2048("g")[1]
        # ya, xa = fu.sky_to_pix(header, ra, dec)
        xa, ya = fu.sky_to_pix(header, ra, dec)
    else:
        s = (to/2 - from_/2)
        yo = s + y
        xo = s + x
        return yo, xo
    print((ya,xa), (yo,xo))
    # return ya, xa
    return xa, ya

            # ra, dec = fu.pix_to_sky(header, self.x, dims-self.y))
    # return y,x

def finding_chart_redo(lens, star_location, slits, savedir=".", runid="123345", typename="lens"):
    figsize((24, 16))
    fontsize(20)
    plt.figure()
    print("Got lens: " + lens.name)
    # ii = to_png([lens.data_512('g'), lens.data_512('r'), lens.data_512('i')], None, band="g", observations=observations)
    # if ii is None:
    # ii = to_png([lens.data_2048('g'), lens.data_2048('r'), lens.data_2048('i')], None, band="g", observations=observations,
    # ii = to_png([lens.data_1024('g'), lens.data_1024('r'), lens.data_1024('i')], None, band="g", observations=observations,
    #             size=1024)
    ii = to_png_redo(lens, star_location, slits)

    plt.imshow(ii)
    scalebar =  ScaleBar(0.263, "arcsec", asec, fixed_value=30, color='blue')

    arrowy = int(ii.size[0] * 0.75)
    arrowx = arrowy
    plt.arrow(arrowy, arrowx, 0, -20, color="green", width=2)
    plt.arrow(arrowy, arrowx, -30, 0, color="green", width=2)
    plt.text(arrowy-4, arrowx-32, "N", color='green', fontsize=20)
    plt.text(arrowy-50, arrowx+5, "E", color='green', fontsize=20)

    plt.gca().add_artist(scalebar)
    plt.text(0, 50, "Run: " + str(runid) + "\nPI: S.Lopez\nTarget=%s (one)\n$\\lambda = 4000-5500\\AA$"%lens.name,
             fontsize=20, color="red")
    plt.title(lens.name)
    # plt.show()
    if savedir is not None:
        plt.savefig(savedir + "/" + lens.name + "_" + typename + "_finding.png", bbox_inches="tight")
    return(ii)

def finding_chart(observations, savedir=".", runid="123345", typename="lens"):
    figsize((24, 16))
    fontsize(20)
    plt.figure()
    lens = observations[0].lens
    print("Got lens: " + lens.name)
    # ii = to_png([lens.data_512('g'), lens.data_512('r'), lens.data_512('i')], None, band="g", observations=observations)
    # if ii is None:
    # ii = to_png([lens.data_2048('g'), lens.data_2048('r'), lens.data_2048('i')], None, band="g", observations=observations,
    ii = to_png([lens.data_1024('g'), lens.data_1024('r'), lens.data_1024('i')], None, band="g", observations=observations,
                size=1024)
    plt.imshow(ii)
    scalebar =  ScaleBar(0.263, "arcsec", asec, fixed_value=30, color='blue')

    arrowy = int(ii.size[0] * 0.75)
    arrowx = arrowy
    plt.arrow(arrowy, arrowx, 0, -20, color="green", width=2)
    plt.arrow(arrowy, arrowx, -30, 0, color="green", width=2)
    plt.text(arrowy-4, arrowx-32, "N", color='green', fontsize=20)
    plt.text(arrowy-48, arrowx+3, "E", color='green', fontsize=20)

    plt.gca().add_artist(scalebar)
    plt.text(0, 50, "Run: " + str(runid) + "\nPI: S.Lopez\nTarget=%s (one)\n$\\lambda = 4000-5500\\AA$"%lens.name,
             fontsize=20, color="red")
    plt.title(lens.name)
    # plt.show()
    if savedir is not None:
        plt.savefig(savedir + "/" + lens.name + "_" + typename + "_finding.png", bbox_inches="tight")
    return(ii)

def hms(ra, dec):
    rahms = ct.decimal_to_hours(ra, dec, return_tuple=False)[0].replace(" ", ":")
    dechms = ct.decimal_to_hours(ra, dec, return_tuple=False)[1].replace(" ", ":")
    return rahms, dechms


def log(string):
    with open("doslit.log", "a") as f:
        f.write(string + "\n")
    print(string)

def fontsize(newsize=None):
    if newsize is not None:
        oldsize = rcParams['font.size']
        if newsize == "bigger":
            newsize = int(oldsize * 1.2)
        elif newsize == "smaller":
            newsize = int(oldsize * 1.2)
        # elif type(newsize) == int:
            # newsize = (int(oldsize * newsize))
        mpl.rcParams['font.size'] = newsize
    return mpl.rcParams['font.size']


def figsize(newsize=None):
    if newsize is not None:
        oldsize = rcParams['figure.figsize']
        if newsize == "bigger":
            newsize = (int(oldsize[0] * 1.2), int(oldsize[1] * 1.2))
        elif newsize == "smaller":
            newsize = (int(oldsize[0] * 1.2), int(oldsize[1] * 1.2))
        elif type(newsize) == int:
            newsize = (int(oldsize[0] * newsize), int(oldsize[1] * newsize))
        mpl.rcParams['figure.figsize'] = newsize
    return mpl.rcParams['figure.figsize']


# TODO
# Cleanup
# Why if fixed to middle is it not finding optimal spots?
# Lens flux threshold?
# OO
# Easier get-flux method
class Lens(object):
    def __init__(self, objname, data, source_mask, lens_mask, header,
            avoid_mask=None, basedir=None):
        self.name = objname
        self.data = data
        self.source_mask = source_mask
        self.lens_mask = lens_mask
        self.avoid_mask = avoid_mask
        self.header = header
        self.basedir = basedir

    def data_512(self, band, flip=True):
        filename = self.basedir + '/fits/512/' + self.name + "_" + band + ".fits"
        return self.data_big(band, flip=flip, filename=filename)

    def data_1024(self, band, flip=True):
        filename = self.basedir + '/1024/' + self.name + "_" + band + ".fits"
        return self.data_big(band, flip=flip, filename=filename)

    def data_2048(self, band, flip=True):
        filename = '/Data/2048/' + self.name + "_" + band + ".fits"
        return self.data_big(band, flip=flip, filename=filename)

    def data_big(self, band, flip=True, filename=None):
        hdulist = fits.open(filename)
        data = hdulist[0].data
        if flip:
            data = np.flipud(data)
        return data, hdulist[0].header


class Template(object):
    def __init__(self, filename):
        self.lines = ct.fal(filename)

    def replace(self, key, value):
        for i, line in enumerate(self.lines):
            if '"' + key + '"' in line:
                self.lines[i] = line.replace('"' + key + '"', 
                        '"' + str(value) + '"')
                return

    def __repr__(self):
        return "\n".join(self.lines)

    def save(self, filename):
        ct.strf(filename, "\n".join(self.lines))


class Position(object):
    def __init__(self, p, theta, radec=None, dims=1024):
        self.x = p[0]
        self.y = p[1]
        self.theta = theta
        self.ra = None
        self.dec = None
        if radec is not None:
            self.ra = radec[0]
            self.dec = radec[1]
        self.dims = dims

    def center(self):
        return self.x, self.y

    def set_ra_dec(self, header):
        self.ra, self.dec = self.get_ra_dec(header)

    def get_ra_dec(self, header):
        if self.ra is not None:
            return (self.ra, self.dec)
        if flipit:
            ra, dec = fu.pix_to_sky(header, self.x, dims-self.y)
        else:
            ra, dec = fu.pix_to_sky(header, self.x, self.y)
        return (ra, dec)

    def from_ra_dec(self, header, ra, dec):
        y, x = fu.sky_to_pix(header, ra, dec)
        if flipit:
            y = self.dims - y
        self.x = x
        self.y = y

    def set_from_ra_dec(self, header):
        self.from_ra_dec(header, self.ra, self.dec)

    def __repr__(self):
        return "(%d, %d), %.2f (%s, %s)" % (self.x, self.y, self.theta, self.ra, self.dec)


def get_slits_from_region_file(filename):
    lines = ct.fal(filename)
    slits = []
    for line in lines:
        if line.startswith("box"):
            toks = line[4:44].split(",")
            # slits.append( Slit(Position(
                # )))
# box(36.90546,-10.65246,1.200",11.000",33.000) # text={XS Slit}
# box(36.90688,-10.65031,1.200",11.000",33.000) # text={XS Slit} color=red

class Slit(object):
    def __init__(self, position, slit_mask):
        self.position = position
        self.slit_mask = slit_mask

    def __repr__(self):
        return ("Slit at %s" % self.position.__repr__())

class Observation(object):
    def __init__(self, lens, best_slit, nod_slits, region_string, star_location=None, times=None):
        self.lens = lens
        self.best_slits = best_slit
        self.nod_slits = nod_slits
        self.region_string = region_string
        self.star_location = star_location
        self.times = times

    def __str__(self):
        return self.lens.name + ": " + str(self.best_slits) + \
               "," + str(self.nod_slits) + ", " + str(self.region_string) \
               + ", " + str(self.star_location) + ", " + str(self.times)
    # return best_slit, nod_slits, output_string
def save_observation(observation, template_name):
    pass


def get_source_flux(lens, slit, bestband="g"):
    slit_source = (slit.slit_mask * lens.source_mask) / 255.0
    sflux = np.sum(lens.data[bestband] * slit_source)
    return sflux


def get_lens_flux(lens, slit, band="r"):
    slit_lens = (slit.slit_mask * lens.lens_mask) / 255.0
    lflux = np.sum(lens.data[band] * slit_lens)
    return lflux

def get_avoid_flux(lens, slit, band="r"):
    slit_lens = (slit.slit_mask * lens.avoid_mask) / 255.0
    lflux = np.sum(lens.data[band] * slit_lens)
    return lflux

def get_slit_mags(lens, slit):
    mags = {"lens": {}, "source": {}}
    for band in "griz":
        mags['lens'][
            band] = 30 - 2.5 * np.log10(get_lens_flux(lens, slit, band=band))
        mags['source'][band] = 30 - 2.5 * np.log10(
            get_source_flux(lens, slit, bestband=band))
        # More than one source?
    return mags


def get_oned_source_slit_fluxes(lens, slit, band):
    slit_source = (slit.slit_mask * lens.source_mask) / 255.0
    # sflux = np.sum(lens.data[bestband] * slit_source)
    flux = lens.data[band] * slit_source
    # print("Flux: ", np.sum(flux), 30 - 2.5 * np.log10(np.sum(flux)))
    im = Image.fromarray(flux)
    im = im.rotate(slit.position.theta)

    total = 0
    sources = []
    data = np.array(im)
    current_source = None
    for x in range(data.shape[0]):
        s = np.sum(data[x, :])
        if s != 0:
            total += s
            if current_source is None:
                current_source = s
            else:
                current_source += s
        else:
            if current_source is not None:
                sources.append(current_source)
                # print("Source:", current_source)
                current_source = None
    # print(total, 30 - 2.5 * np.log10(total))
    return np.array(sources), flux, np.array(im)


def embed(image, new_dim):
    bigger = np.zeros((new_dim, new_dim))
    smallsize = image.shape[0]
    h = new_dim / 2
    s = smallsize / 2
    bigger[h - s:h + s, h - s:h + s] = image
    return bigger


def distance(ax, ay, bx, by):
    return math.sqrt((by - ay)**2 + (bx - ax)**2)


def rotated_about(ax, ay, bx, by, angle):
    radius = distance(ax, ay, bx, by)
    angle += math.atan2(ay - by, ax - bx)
    return (round(bx + radius * math.cos(angle)),
            round(by + radius * math.sin(angle)))

def oned_flux(lens, slit, band="g", mask=None):
    slit_ap = slit.slit_mask
    flux = lens.data[band] * slit_ap
    if mask is not None:
        smask = (slit_ap * mask)/255
        flux = lens.data[band] * smask
    if flux.sum() == 0:
#         return np.zeros(np.array(slit.slit_mask.shape).max())
        return np.zeros(int(11/.263))
    im = Image.fromarray(flux)
    im = im.rotate(slit.position.theta)
    ap = Image.fromarray(255.0*slit_ap)
    ap = ap.rotate(slit.position.theta)

    nonzero = np.where(np.array(ap) != 0)
    mins = min(nonzero[0]), min(nonzero[1])
    maxes = max(nonzero[0]), min(nonzero[1])
    slit_only = np.array(im)[mins[0]:maxes[0], :]
    od_flux = slit_only.sum(axis=1)
    return od_flux

def old_oned_flux(lens, slit, band="g"):
    slit_ap = slit.slit_mask
    flux = lens.data[band] * slit_ap
    im = Image.fromarray(flux)
    im = im.rotate(slit.position.theta)
    nonzero = np.where(np.array(im) != 0)
    mins = min(nonzero[0]), min(nonzero[1])
    maxes = max(nonzero[0]), min(nonzero[1])
    slit_only = np.array(im)[mins[0]:maxes[0], :]
    oned_flux = slit_only.sum(axis=1)
    return oned_flux


def rotate(degrees, sw, sh, dims, rect_center=None, image=None, fill=255, outline=None):
    if image is None:
        image = Image.new('1', (dims, dims))
    draw = ImageDraw.Draw(image)

    if rect_center is None:
        rect_center = (dims / 2, dims / 2)

    rect_vertices = ((rect_center[0] + sw / 2, rect_center[1] + sh / 2),
                     (rect_center[0] + sw / 2, rect_center[1] - sh / 2),
                     (rect_center[0] - sw / 2,
                      rect_center[1] - sh / 2), (rect_center[0] - sw / 2,
                                                 rect_center[1] + sh / 2))

    rect_vertices = [
        rotated_about(x, y, rect_center[0], rect_center[1],
                      math.radians(degrees)) for x, y in rect_vertices
    ]

    draw.polygon(rect_vertices, fill=fill, outline=outline)
    return image


def region_file(image, center, sh, sw, angle):
    pass


def get_lens(objname, maskdir=".", fitsdir="."):
    mask_lens = np.array(Image.open(maskdir+ "/masks_lens/" + objname + ".bmp"))
    mask_source = np.array(
        Image.open(maskdir + "/masks_source/" + objname + ".bmp"))
    avoid_mask = None
    if os.path.exists(maskdir+ "/masks_avoid/" + objname + ".bmp"):
        avoid_mask = np.array(Image.open(maskdir + "/masks_a/" + objname + ".bmp"))
    if not flipit:
        mask_lens = np.flipud(mask_lens)
        mask_source = np.flipud(mask_source)
        if avoid_mask is not None:
            avoid_mask = np.flipud(avoid_mask)

    bands = "griz"
    data = {}
    header = None
    for band in bands:
        filename = fitsdir + '/100/' + objname + "_" + band + ".fits"
        # print(filename)
        if os.path.exists(filename):
            hdulist = fits.open(filename)
            data[band] = hdulist[0].data
            if flipit:
                data[band] = np.flipud(data[band])
                # This is an artifact of png creation; they are flipped, masks
                # are made from pngs
            if band == "g":
                header = hdulist[0].header

    lens = Lens(objname, data, mask_source, mask_lens, header, avoid_mask=avoid_mask, basedir=fitsdir)
    return lens

def process_star(lens, star, ra, dec):
    closest_ra, closest_dec = star['ra'], star['dec']
    regions = ""

    box_size = 19
    d = lens.data_1024('r', flip=False)
    y, x = fu.sky_to_pix(d[1], closest_ra, closest_dec)
    b = int(box_size / 2)
    star = d[0][x - b:x + b, y - b:y + b]
    if star.shape[0] != box_size - 1 or star.shape[1] != box_size - 1:
        print("Star off edge, not recentering.")
        center_ra, center_dec = closest_ra, closest_dec
    else:
        if False:
            center = np.unravel_index(star.argmax(), star.shape)
        else:
            center = centroid(star)
        center_ra, center_dec = fu.pix_to_sky(d[1], center[1] + (y - b),
                                              center[0] + (x - b))
    regions += "line(%.6f,%.6f,%.6f,%.6f) # line=0 0\n" % (center_ra,
                                                           center_dec, ra, dec)
    return regions, (center_ra, center_dec)

def optimise_observation(lens, # dictionary
                    lens_flux_minimum=1000,
                    bestband="g",
                    nod_flux=0.95,
                    make_diag=False,
                    centered=False,
                    force_angle=None,
                    target="lens",
                    # lens_only=False,
                    # source_only=False,
                    avoid=0,
                    dataframe=None,
                    max_distance=30,
                    times=None,
                    slitcolors=["cyan", "magenta"],
                    source_flux_maximum=0,
                    do_star=True,
                         just_slits=False,
                         clear_stars_cache=False,
                         **kwargs):
    # Open the fits files
    # Open the masks
    # rotate slit to optimise source flux
    center = (dims / 2, dims / 2)
    # if trim_slit:
    # slitheight[0] = 34

    source_flux = 0
    lens_flux = 0
    target_flux = 0
    best_position = None
    best_slit = None

    iteration = 0
    # if centered:
    #     x_range = [50]
    #     y_range = [50]
    # else:

    x_range = range(50 - max_distance, 50 + max_distance, 2)
    y_range = range(50 - max_distance, 50 + max_distance, 2)

    theta_range = range(0, 180, 5)

    if force_angle is not None:
        theta_range = [180 - force_angle]
    its = 0
    for x in x_range:
        for y in y_range:
            for theta in theta_range:
                its += 1
                # ct.progbar(iteration, ((dims - 60) / 2)**2)
                center = (x, y)

                slit = Slit(
                    Position(center, theta),
                    np.array(
                        rotate(
                            theta,
                            slitwidth,
                            slitheight[0],
                            dims,
                            rect_center=center)))
                # Comes back as boolean - is pixel in slit?
                if target == "lens":
                    flux = get_lens_flux(lens, slit)
                else:
                    flux = get_source_flux(lens, slit, bestband=bestband)


                if avoid > 0:
                    aflux = get_avoid_flux(lens, slit)

                update = False
                if False: # before
                    if lens_only:
                        if sflux <= source_flux_maximum and lflux > lens_flux:
                            update = True

                    elif (sflux > source_flux) and (lflux >= lens_flux_minimum):
                        if avoid == 0 or aflux <= avoid:
                            update = True
                            # print("<- best")
                if flux > target_flux:
                    update = True

                if update:
                    target_flux = flux
                    # source_flux = sflux
                    # lens_flux = lflux
                    best_position = Position(center, theta)
                    best_slit = slit

            iteration += 1

    print("Iterations of loop:", its)
    # x, y = best_position[0]
    print(best_position)
    if best_position is None:
        return None
    if just_slits:
        return best_position

    iterations = 0
    # mask_lens, mask_source)
    # print(positions)

    # Convert center to ra,dec
    rf = ""  # region file
    # rf = make_region_file(best_position, header)
    nod_slits = do_nod2(
        lens, best_slit, minimum_retained=nod_flux, save_all=False, avoid=avoid, max_dist=11)
    for ss in nod_slits:
        # print("**** NOD POS:", ss.position)
        ss.position.set_ra_dec(lens.header)
        # print("**** NOD POS:", ss.position)

    rf += "\n" + make_region_file(nod_slits[0].position, lens.header, color=slitcolors[0], comment="XS Slit A")
    rf += "\n" + make_region_file(
        nod_slits[1].position, lens.header, color=slitcolors[1], comment="XS Slit B")
    # ra, dec = fu.pix_to_sky(lens.header, best_slit.position.x, best_slit.position.y)
    if flipit:
        ra, dec = fu.pix_to_sky(lens.header, nod_slits[0].position.x,
                                dims - nod_slits[0].position.y)
        ra2, dec2 = fu.pix_to_sky(lens.header, nod_slits[1].position.x,
                                  dims - nod_slits[1].position.y)
    else:
        ra, dec = fu.pix_to_sky(lens.header, nod_slits[0].position.x,
                                nod_slits[0].position.y)
        ra2, dec2 = fu.pix_to_sky(lens.header, nod_slits[1].position.x,
                                  nod_slits[1].position.y)
    if clear_stars_cache:
        try:
            os.remove(f"stars/{lens.name}.csv")
        except:
            pass

    try:
        stars = pd.read_csv(f"stars/{lens.name}.csv")
    except Exception as e:
        mag_min = 18
        mag_max = 20
        distance = 1

        if 'star_min_mag' in kwargs:
            mag_min = kwargs['star_min_mag']
        if "star_max_mag" in kwargs:
            mag_max = kwargs['star_max_mag']
        if "star_max_dist" in kwargs:
            distance = kwargs['star_max_dist']
        skip_stars = kwargs['skip_stars']

        stars = starutils.do_stars(lens, ra, dec, distance=distance, mag_min=mag_min, mag_max=mag_max)
        if stars is not None and len(stars) <= skip_stars:
            stars = None

        if stars is not None:
            stars = stars[skip_stars:]
            stars.to_csv(f"stars/{lens.name}.csv")
  
    if do_star and stars is not None: 
        closest, cdistance = starutils.get_closest(stars, ra, dec)
        stars = process_star(lens, stars.loc[closest], ra, dec)
    else:
        do_star = False
        if stars is None:
            log(f"No star for {lens.name}")

    if do_star:
        rf += "\n" + stars[0]
        output_string = "%s,%.6f,%.6f,%.6f,%.6f,%.2f\n" % (lens.name, ra, dec, stars[1][0],
                                                       stars[1][1], 180 - nod_slits[0].position.theta)
    else:
        output_string = "%s,%.6f,%.6f,%s,%s,%.2f\n" % (lens.name, ra, dec, "-",
                                                       "-", 180 - nod_slits[0].position.theta)

    if dataframe is not None:

        dataframe.loc[lens.name, "ra"] = ra
        dataframe.loc[lens.name, "dec"] = dec
        dataframe.loc[lens.name, "ra_hms"] = hms(ra, dec)[0]
        dataframe.loc[lens.name, "dec_dms"] = hms(ra, dec)[1]
        dataframe.loc[lens.name, "star_ra"] = stars[1][0]
        dataframe.loc[lens.name, "star_dec"] = stars[1][1]
        dataframe.loc[lens.name, "star_ra_hms"] = hms(stars[1][0], stars[1][1])[0]
        dataframe.loc[lens.name, "star_dec_dms"] = hms(stars[1][0], stars[1][1])[1]
        # dataframe.loc[lens.name, "star_ra"] = stars[1][0]
        # dataframe.loc[lens.name, "star_dec"] = stars[1][1]
        dataframe.loc[lens.name, "offset_ra"] = -3600 * ((stars[1][0] - ra) * np.cos(stars[1][1] * np.pi / 180.))
        dataframe.loc[lens.name, "offset_dec"] = -3600 * (stars[1][1] - dec)
        pa = 180 - nod_slits[0].position.theta
        dataframe.loc[lens.name, "pa"] = pa
        dataframe.loc[lens.name, "nod_distance"] = 3600 * np.sqrt((dec2 - dec) ** 2 + \
                                                                  (np.cos(dec * np.pi / 180.) * (ra2 - ra)) ** 2)
        if (pa <= 90 and (dec2 - dec) < 0) or (pa > 90 and (dec2 - dec) > 0):
            dataframe.loc[lens.name, "nod_distance"] *= -1

    # save_region_file(basedir + '/' + lens.name + ".reg", rf)
    if make_diag:
        make_diagnostic(lens, nod_slits)
    if do_star:
        obs = Observation(lens, best_slit, nod_slits, output_string,
                      star_location=stars[1], times=times)
    else:
        obs = Observation(lens, best_slit, nod_slits, output_string,
                      None, times=times)

    pkl.dump(obs, open(lens.name + ".pkl", "wb"))
    return dict(slit=best_slit, nod=nod_slits, string=output_string, obs=obs, rf=rf)


def center_star(data, x, y, box_size=7):
    b = int(box_size / 2)
    star = data[x - b:x + b, y - b:y + b]
    center = np.unravel_index(star.argmax(), star.shape)
    return center


def distance_sky_px(start, finish, header):
    ra1, dec1 = fu.pix_to_sky(header, start[0], start[1])
    ra2, dec2 = fu.pix_to_sky(header, finish[0], finish[1])
    return distance_sky((ra1, dec1), (ra2, dec2))


def distance_sky(start, finish):
    delta_ra = np.abs(np.cos(finish[1] * np.pi / 180) * (start[0] - finish[0]))
    delta_dec = np.abs(start[1] - finish[1])
    distance = np.sqrt(delta_ra**2 + delta_dec**2)
    return distance
    # return delta_ra, delta_dec


def offset_sky(start, finish, arcsec=True):
    delta_ra = np.abs(np.cos(finish[1] * np.pi / 180) * (start[0] - finish[0]))
    delta_dec = np.abs(start[1] - finish[1])
    if arcsec:
        return delta_ra * -3600, delta_dec * -3600
    else:
        return delta_ra, delta_dec


def get_offset(objid):
    lines = ct.fal(basedir + "/" + objid + "_stars.reg")
    for line in lines:
        if ("line") in line:
            # line(22.51697,-37.72629,22.51181,-37.74937)
            ra = float(line.split(",")[0].split("(")[1])
            dec = float(line.split(",")[1])
            ra2 = float(line.split(",")[2])
            dec2 = float(line.split(",")[3].split(")")[0])
            return offset_sky((ra, dec), (ra2, dec2))


def save_region_file(filename, regionfile):
    regionfile = \
"""# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
""" + regionfile
    ct.strf(filename, regionfile)


def make_region_file(best_position, header, comment="XS Slit", color=None):
    # Flipud
    if flipit:
        best_position = Position((best_position.x, dims - best_position.y),
                                 best_position.theta)
    ra, dec = fu.pix_to_sky(header, best_position.x, best_position.y)
    regionfile = ""
    colorstring = ""
    if color is not None:
        colorstring = " color=%s" % color
    regionfile += "box(%.6f,%.6f,1.200\",11.000\",%.3f) # text={%s}%s" % (
        ra, dec, 180 - best_position.theta, comment, colorstring)
    # print(regionfile)
    return regionfile


def dy(distance, m):
    return m * dx(distance, m)


def dx(distance, m):
    return np.sqrt(distance / (m**2 + 1))


def do_nod2(lens, slit, save_all=False, minimum_retained=.95, avoid=0, max_dist=9, min_dist=2):
    slits = [slit, slit]
    lflux_target = get_lens_flux(lens, slit)
    sflux_target = get_source_flux(lens, slit)

    print("Nodding: ", slit.position, lflux_target, sflux_target)

    iterations = 0
    startx, starty = slit.position.center()
    x, y = startx, starty
    alpha = (np.pi / 180) * (90 - slit.position.theta)
    m = np.tan(alpha)
    # print(
    # "The line is: y = %.2fx + %.2f" % (m, starty - startx * np.tan(alpha)))

    direction = 0
    while True:
        # iterations < 100 and
        distance = 3600*distance_sky_px((startx, startx), (x, y), lens.header)
        # distance_px = np.sqrt((x-startx)**2 + (y-starty)**2)
        # distance_compare = distance_px * .263
        # print("Distance: ", distance, distance_px, distance_compare)
        if distance > max_dist or iterations > 100:
            if direction == 1:
                break
            else:
                direction = 1
                iterations = 0
                x, y = startx, starty
                continue

        if direction == 0:
            new_x, new_y = (x + dx(1, m), y - dy(1, m))
        else:
            new_x, new_y = (x - dx(1, m), y + dy(1, m))

        if distance < min_dist:
            x, y = new_x, new_y
            iterations += 1
            continue

        new_slit = Slit(
            Position((new_x, new_y), slit.position.theta),
            np.array(
                rotate(
                    slit.position.theta,
                    slitwidth,
                    slitheight[0],
                    dims,
                    rect_center=(new_x, new_y))))
        lflux = get_lens_flux(lens, new_slit)
        sflux = get_source_flux(lens, new_slit)
        aflux = 0
        if avoid > 0:
            aflux = get_avoid_flux(lens, new_slit)

        # print("Working:", sflux, sflux_target)
        if sflux >= sflux_target * minimum_retained and lflux >= lflux_target * minimum_retained:
            if aflux <= avoid:
                slits[direction] = new_slit
        if save_all:
            slits.append(new_slit)
        x, y = new_x, new_y
        iterations += 1

    return slits


# def do_nod(slit, position, data, objid, header, threshold=100):
# x, y = position[0]
# angle = position[1]
# arcsec = 1 / .263
# distance = 2 * arcsec
# flux = 1e6
# while flux > threshold and distance < 30:
# dy = distance * np.sin((2 * np.pi / 360) * angle)  # degrees
# dx = distance * np.cos((2 * np.pi / 360) * angle)
# new_x, new_y = x + dx, y + dy
# new_slit = np.array(
# rotate(
# angle, slitwidth, slitheight, dims, rect_center=(new_x,
# new_y)))
# flux = np.sum(data * new_slit)
# x, y = new_x, new_y
# distance += 1
# rf = make_region_file(
# ((x, y), angle),
# header,
# comment="Nod slit %.1f\"" % (distance / arcsec))
# return ((x, y), angle), new_slit, rf

def source_area(lens, slit):
    source_area = np.where(slit.slit_mask*lens.source_mask/255*lens.data['g'] > 0)
    return len(source_area[0])


def make_diagnostic(lens, slits):
    table, label = make_texts(lens, slits[0])
    table += "\nSource area: %d, %d" % (source_area(lens, slits[0]), source_area(lens, slits[1]))
    ct.strf(lens.name + "_slitmags.txt", table)
    figsize((12, 8))
    png = Image.open("/home/coljac/phd/Proposals/xshooter/working/" +
                     lens.name + ".png")
    slit = slits[0]
    slit2 = slits[1]
    plt.figure()
    plt.title(lens.name)
    plt.subplot(1, 3, 1)
    plt.imshow(png)
    plt.subplot(1, 3, 2)
    plt.imshow(slit.slit_mask * 255 + np.array(lens.source_mask) +
               np.array(lens.lens_mask) * 100)
    plt.subplot(1, 3, 3)
    plt.imshow(slit2.slit_mask * 255 + np.array(lens.source_mask) +
               np.array(lens.lens_mask) * 100)
    plt.text(-240, 180, label, fontsize=20)
    plt.savefig(lens.name + "_slits.png", bbox_inches="tight")
    plt.figure()
    # colors = ["red", "blue"]
    a = np.zeros(47)
    b = np.zeros(47)


    slit = slits[0]
    slit2 = slits[1]
    for band in "gr":
        plt.figure()
        oned = oned_flux(lens, slit, band=band, mask=lens.source_mask)
        plt.plot(np.arange(oned.shape[0]), oned, "b-", label="src")
        oned = oned_flux(lens, slit, band=band, mask=lens.lens_mask)
        plt.plot(np.arange(oned.shape[0]), oned, "r-", label="lens")

        oned = oned_flux(lens, slit2, band=band, mask=lens.source_mask)
        plt.plot(np.arange(oned.shape[0]), oned, "b--", label="src nod")
        oned = oned_flux(lens, slit2, band=band, mask=lens.lens_mask)
        plt.plot(np.arange(oned.shape[0]), oned, "r--", label="lens nod") 
        plt.title(lens.name + " in " + band)
        plt.legend()
        plt.savefig(lens.name + "_slitflux_" + band + ".png", bbox_inches="tight")

def old_slitflux():
    for band in "gri":
        oned = oned_flux(lens, slit, band=band)
        a[0:oned.shape[0]] = oned
        oned2 = oned_flux(lens, slit2, band=band)
        b[0:oned2.shape[0]] = oned2
        plt.plot(np.arange(a.shape[0]), a - b)
    plt.savefig(lens.name + "_slitflux.png", bbox_inches="tight")


def make_texts(lens, slit):
    # Slit mags
    t = PrettyTable()
    t.field_names = ["band", "lens", "source"]
    mags = get_slit_mags(lens, slit)
    fig_label = "Lens: "
    for band in "griz":
        t.add_row(
            [band,
             "%.2f" % mags['lens'][band],
             "%.2f" % mags['source'][band]])
        fig_label += "$m_%s = %.1f$, " % (band, mags['lens'][band])
    fig_label += "\nSources: "
    text_output = "Slit mags:\n" + t.get_string() + "\n"
    t = None
    smags = None
    j = 0
    for band in "griz":
        sf = get_oned_source_slit_fluxes(lens, slit, band)[0]
        if t is None:
            smags = np.zeros((4, len(sf)))
            t = PrettyTable()
            t.field_names = ["band"] + [
                'source ' + str(i + 1) for i in range(len(sf))
            ]
        t.add_row([band] + ["%.2f" % (30 - 2.5*np.log10(sf[i])) for i in range(sf.shape[0])])
        for i in range(sf.shape[0]):
            smags[j, i] = 30-2.5*np.log10(sf[i])
        j += 1
    for i in range(smags.shape[1]):
        for j in range(4):
            fig_label += "$m^%d_%s = %.1f$, " % (i+1, "griz"[j], smags[j, i])
        fig_label += "\n"
    text_output += "\nIndividual sources:\n" + t.get_string()

    # Exposure times
    t = PrettyTable()
    t.field_names = ["", "time"]

    return text_output, fig_label


def centroid(data):
    # x, y = centroid_com(data)
    # return x,y
    x1, y1 = centroid_2dg(data)
    # return int(np.rint(y1)), int(np.rint(x1))
    return y1, x1


def make_obx(o1, o2, output, summary_file):
    summary = pd.read_csv(summary_file, index_col="name")
    row = summary.loc[o1.lens.name]
    if o2 is None:
        template = Template("single.obx")
    else:
        template = Template("double.obx")
    nod_throw = row.nod_distance
    if row.nod_distance < 0:
        o1.nod_slits = [o1.nod_slits[1], o1.nod_slits[0]]
        nod_throw *= -1

    # if o2 is None:
        # coordsA = o1.nod_slits[0].position.get_ra_dec(o1.lens.header)
        # coordsB = o1.nod_slits[1].position.get_ra_dec(o1.lens.header)
    # else:
    coordsA = o1.nod_slits[0].position.get_ra_dec(o1.lens.header)
    coordsB = o1.nod_slits[1].position.get_ra_dec(o1.lens.header)

    coordsAB = (coordsA[0] + coordsB[0])/2, (coordsA[1]+coordsB[1])/2
    offsetAB = calc_offset(o1.star_location, coordsAB)

    if True:
        star_ra_dec = hms(*o1.star_location)
        template.replace("TARGET-NAME", o1.lens.name)
        # template.replace("RA", "%s" % row.ra_hms)
        # template.replace("DEC", "%s" % row.dec_dms)
        template.replace("OBSERVATION-NAME", "Lens and source in one")
        template.replace("FINDING-CHART", o1.lens.name + "_finding.jpg")
        template.replace("STAR-RA", "%s" % row.star_ra_hms)
        template.replace("STAR-DEC", "%s" % row.star_dec_dms)
        template.replace("POSITION-ANGLE", "%s" % row.pa)

        template.replace("OFFSETALPHA", "%.3f" % offsetAB[0])
        template.replace("OFFSETDELTA", "%.3f" % offsetAB[1])
        # if o2 is None:
            # template.replace("OFFSETALPHA", "%.3f" % row.offset_ra)
            # template.replace("OFFSETDELTA", "%.3f" % row.offset_dec)
        template.replace("NOD THROW", "%.2f" % nod_throw)
    if o2 is None:
        template.replace("OBS-NAME", o1.lens.name + " single obs")
    if o2 is not None: 

        template.replace("OBS-NAME", o1.lens.name + " double obs")
        template.replace("OFFSETALPHA", "%.3f" % offsetAB[0])
        template.replace("OFFSETDELTA", "%.3f" % offsetAB[1])

        template.replace("OBS-NAME", o1.lens.name + " double obs")
        coords1 = o1.nod_slits[0].position.get_ra_dec(o1.lens.header)
        coords2 = o2.nod_slits[0].position.get_ra_dec(o1.lens.header)
        coords3 = o2.nod_slits[1].position.get_ra_dec(o1.lens.header)
        offsetAB_to_A = calc_offset(coordsAB, coords2)
        # offset = calc_offset(coords1, coords2)
        offset2 = calc_offset(coords2, coords3)
        template.replace("RELOFF1", "%.2f %.2f %.2f %.2f" % (offsetAB_to_A[0], offset2[0], -1*offset2[0], offset2[0]))
        template.replace("RELOFF2", "%.2f %.2f %.2f %.2f" % (offsetAB_to_A[1], offset2[1], -1*offset2[1], offset2[1]))
    reg_patch = "point(%.6fd,%.6fd) # color=blue\n" % (coordsAB[0], coordsAB[1])
    reg_patch += "line(%.6fd,%.6fd,%.6f,%.6f) # line=0 0\n" % (o1.star_location[0], 
                o1.star_location[1], coordsAB[0], coordsAB[1])
    mpname = "mp.reg" if o2 is not None else "mp1.reg"
    ct.strf(o1.lens.name + mpname, reg_patch)

    template.save(output)


def calc_offset(coords1, coords2):
    offset_ra = 3600 * (coords2[0] - coords1[0]) * np.cos(coords2[1]*np.pi/180)
    offset_dec = 3600 * (coords2[1] - coords1[1])
    print(offset_dec, coords2[1],"-", coords1[1])
    return offset_ra, offset_dec


def obs(argv, make_finding):
    ids = ct.fal(argv[1])
    dir1 = argv[2]
    dir2 = None
    if len(argv) > 3:
        dir2 = argv[3]
    o2 = None
    count = 0
    for i in ids:
        count += 1
        print("Making obs for ", i, "%d/%d" % (count, len(ids)))
        o1 = pkl.load(open(dir1 + "/" + i + ".pkl", "rb"))
        if dir2 is not None:
            o2 = pkl.load(open(dir2 + "/" + i + ".pkl", "rb"))
            print(i, o1.nod_slits[0].position.get_ra_dec(o1.lens.header))
        # make_obx(o1, o2, (dir2 if dir2 is not None else dir1) + "/" + i + ".obx", dir1 + "/summary.csv")
        make_obx(o1, o2, dir1 + "/" + i + ".obx", dir1 + "/summary.csv")
        if make_finding:
            if o2 is None:
                finding_chart([o1], savedir=dir1)
            else:
                finding_chart([o1, o2], savedir=dir2)

def oldmain(argv):
    centered = ct.argb(argv, "-c")
    nod_flux = ct.argv(argv, "-n")
    lens_only = ct.argb(argv, "-l")
    bcolors = ct.argb(argv, "-b")
    # if(bcolors):
        # slitcolors[0] = "cyan"
        # slitcolors[1] = "magenta"
    source_flux_maximum = 0
    if nod_flux is None:
        nod_flux = 0.90 # 0.9
    else:
        nod_flux = float(nod_flux)
    angle = ct.argv(argv, "-a")
    try:
        angle = float(angle)
    except:
        pass
    avoid = ct.argv(argv, "-v")
    if avoid is None:
        avoid = 0
    else:
        avoid = float(avoid)
    sfm = ct.argv(argv, "-x")
    if sfm is not None:
        sfm = float(sfm)
    else:
        sfm = source_flux_maximum

    ids = ct.fal(argv[1])
    bb = "g"
    if len(argv) > 2:
        bb = argv[2]
    lfm = 1000
    if len(argv) > 3:
        lfm = int(argv[3])
    # angle = None
    # if len(argv) > 4:
        # angle = float(argv[4])
    # avoid = 0
    # if len(argv) > 5:
        # avoid = float(argv[5])

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    print("Running at %s with bb=%s, lfm=%d, angle=%s, avoid=%s, nod_flux=%.2f, centered=%s, sfm=%.1f" % (now, bb, lfm, str(angle), str(avoid), nod_flux, str(centered), source_flux_maximum))
    log("Running at %s with bb=%s, lfm=%d, angle=%s, avoid=%s, nod_flux=%.2f, centered=%s, sfm=%.1f" % (now, bb, lfm, str(angle), str(avoid), nod_flux, str(centered), source_flux_maximum))
    summary = None
    if os.path.exists("summary.csv"):
        summary = pd.read_csv("summary.csv", index_col="name")
    else:
        summary = pd.DataFrame(columns=["name","ra","dec","ra_hms","dec_dms","star_ra","star_dec","star_ra_hms","star_dec_dms",
            "offset_ra", "offset_dec", "pa", "nod_distance"])
        summary = summary.set_index("name")
    # outputs = "name,ra,dec,star_ra,star_dec,pa\n"

    for i in ids:
        bestband = bb
        if i == "DESJ2249-0110" or i == "DESJ2041-5514":
            bestband = "r"
        if os.path.exists(i + ".reg"):
            # print("Skip")
            continue
        print("Processing: " + str(i))
        log("Target: " + i)
        lens = get_lens(i)
        if not lens_only:
            slitcolors = ["yellow" ,"red"]
            obs_src = optimise_observation(lens, lens_flux_minimum=lfm, nod_flux=nod_flux,
                bestband="g", target="source", centered=centered, force_angle=angle,
                slitcolors=slitcolors,
                avoid=avoid, dataframe=summary, source_flux_maximum=sfm, **kwargs)
            angle = 180 - obs_src['slit'].position.theta
            slitcolors = ["cyan" ,"magenta"]
            obs_lens = optimise_observation(lens, lens_flux_minimum=lfm, nod_flux=nod_flux,
                bestband="r", centered=centered, force_angle=angle,
                slitcolors=slitcolors,
                avoid=avoid, dataframe=summary, do_star=False,
                source_flux_maximum=sfm, **kwargs)
        else:
            obs_lens = optimise_observation(lens, lens_flux_minimum=lfm, nod_flux=nod_flux,
                bestband=bestband, centered=centered, force_angle=angle,
                avoid=avoid, dataframe=summary, lens_only=lens_only,
                source_flux_maximum=sfm, **kwargs)
        obs = obs_lens if lens_only else obs_src

        if obs is None:
            print("Failed with " + i + ", lowering threshold")
            lfm /= 2
            obs = optimise_observation(lens, lens_flux_minimum=lfm, nod_flux=nod_flux, 
                    bestband=bestband, centered=centered, force_angle=angle, avoid=avoid, 
                    lens_only=lens_only, dataframe=summary, source_flux_maximum=sfm, **kwargs)
        # outputs += obs[2]
        if obs is None:
            print("Warning**: Null observation for ", lens.name)
            log("Warning**: Null observation for " + lens.name)
        else:
            print("RETURNED:", type(obs))
            for k, v in obs.items():
                print(f"{k}: {v}")
            save_region_file(f"{i}.reg", obs_lens['rf'] + "\n" +obs_src['rf'])
            # print(obs['obs'].region_string)
            # print(obs[1])
            # print(obs[2])
    # ct.strf("summary.csv", outputs)
    summary.to_csv("summary.csv")
    ct.strf("README.md", "Configuration: Lens flux minimum=%d, best band = %s, centered=%s" % (lfm, bb, str(centered)))



def de(s, name):
    if not os.path.exists(s):
        print(f"{name} {s} not found.")
        raise Exception 
    return True

def do_slits(catalog, position_angle=None, fits_dir=".", masks_dir=".", lens_only=False, source_only=False, nod_flux=.9,
        source_flux_maximum=0, lens_flux_minimum=1000, target=None, no_cache=False, clear_stars_cache=False, **kwargs):

    cat = pd.read_csv(catalog, index_col="name")
    for i, row in cat.iterrows():
        if target is not None and i != target:
            continue
            
        log("Target: " + i)
        if os.path.exists(i + ".reg") and not no_cache:
            print(f"Region exists for {i}")
            log("Region exists.")
            continue
        try:
            lens = get_lens(i, maskdir=masks_dir, fitsdir=fits_dir)
        except Exception as e:
            log(f"Can't do {i}, missing data: {repr(e)}")
            continue


        print("Processing: " + str(i))
        angle = None if position_angle < 0 else position_angle
        region_string = ""
        if not lens_only:
            slitcolors = ["yellow" ,"red"]
            obs_src = optimise_observation(lens, lens_flux_minimum=lens_flux_minimum, nod_flux=nod_flux,
                bestband="g", target="source", 
                # centered=centered, 
                force_angle=angle,
                slitcolors=slitcolors,
                # avoid=avoid, 
                # dataframe=summary, 
                source_flux_maximum=source_flux_maximum, 
                clear_stars_cache=clear_stars_cache,
                **kwargs)
            angle = 180 - obs_src['slit'].position.theta
            region_string += obs_src['rf'] + "\n"

        if not source_only:
            slitcolors = ["cyan" ,"magenta"]
            obs_lens = optimise_observation(lens, lens_flux_minimum=lens_flux_minimum, nod_flux=nod_flux,
                bestband="r", 
                # centered=centered, 
                force_angle=angle,
                slitcolors=slitcolors,
                # avoid=avoid, 
                # dataframe=summary, 
                do_star=False,
                source_flux_maximum=source_flux_maximum,
                clear_stars_cache=clear_stars_cache, **kwargs)
            region_string += obs_lens['rf'] + "\n"        
        save_region_file(f"{i}.reg", region_string)
        # else:
            # obs_lens = optimise_observation(lens, lens_flux_minimum=lfm, nod_flux=nod_flux,
                # bestband=bestband, centered=centered, force_angle=angle,
                # avoid=avoid, dataframe=summary, lens_only=lens_only,
                # source_flux_maximum=sfm)
        # obs = obs_lens if lens_only else obs_src

        # if obs is None:
            # print("Failed with " + i + ", lowering threshold")
            # lfm /= 2
            # obs = optimise_observation(lens, lens_flux_minimum=lfm, nod_flux=nod_flux, 
                    # bestband=bestband, centered=centered, force_angle=angle, avoid=avoid, 
                    # lens_only=lens_only, dataframe=summary, source_flux_maximum=sfm)
        # outputs += obs[2]
        # if obs is None:
            # print("Warning**: Null observation for ", lens.name)
            # log("Warning**: Null observation for " + lens.name)
        # else:
            # print("RETURNED:", type(obs))
            # for k, v in obs.items():
                # print(f"{k}: {v}")
            # save_region_file(f"{i}_test.reg", obs_lens['rf'] + "\n" +obs_src['rf'])angle,
            # print(obs['obs'].region_string)
            # print(obs[1])
            # print(obs[2])
    # ct.strf("summary.csv", outputs)
    # summary.to_csv("summary.csv")
    # ct.strf("README.md", "Configuration: Lens flux minimum=%d, best band = %s, centered=%s" % (lfm, bb, str(centered)))


def main(argv):
    parser = argparse.ArgumentParser(description="XShooter slits for lenses")
    parser.add_argument(dest="catalog", help="Catalogue of targets", nargs="?")
    parser.add_argument('--position-angle', "-p", type=int, default=-1, help="Fix position angle")
    parser.add_argument('--make-finding-charts', "-r", action="store_true", help="Make finding charts only")
    parser.add_argument('--fits-dir', "-f", type=str, required=False, default="./fits")
    parser.add_argument('--masks-dir', "-m", type=str, required=False, default="./masks")
    parser.add_argument('--clear-stars-cache', "-x", action="store_true", help="Don't use cached stars")
    parser.add_argument('--no-cache', "-c", action="store_true", help="Don't use cached optimisations")
    parser.add_argument('--lens-only', "-l", action="store_true")
    parser.add_argument('--source-only', "-s", action="store_true")
    parser.add_argument('--nod-flux', type=float, default=0.9)
    parser.add_argument('--source-flux-maximum', type=float, default=0.0)
    parser.add_argument('--lens-flux-minimum', type=float, default=1000)
    parser.add_argument('--target', '-t', type=str, required=False, help="Target name to run (still uses catalogue)")

    parser.add_argument('--star-min-mag', type=float, default=18.0, help="Minimum star magnitude")
    parser.add_argument('--star-max-mag', type=float, default=20.0, help="Maximum star magnitude")
    parser.add_argument('--star-max-dist', type=float, default=1.0, help="Maximum star distance (arcmin)")
    parser.add_argument('--skip-stars', type=int, default=0, help="In case first star no good")

    args = parser.parse_args()
    with open("doslit.log", "w") as f:
        pass

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    log(f"Running at {now}, with settings: {vars(args)}")
    
    try:
        de(args.fits_dir, "Fits dir") 
        de(args.masks_dir, "Mask dir") 
        de(args.catalog, "Catalogue") 
    except:
        return

    do_slits(**(vars(args)))

    
    

if __name__ == "__main__":
    main(sys.argv)

