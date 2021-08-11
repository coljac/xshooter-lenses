import os
import sys
import numpy as np
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
# import fitsutils as fu
from astroquery.vizier import Vizier
from astroquery.vo_conesearch import ConeSearch
from astroquery.vo_conesearch import conesearch
from astroquery.gaia import Gaia
from astroquery import esasky
from prettytable import PrettyTable
from astroquery.mast import Observations
from astropy.coordinates import SkyCoord
# from icecream import ic
import coltools as ct
import warnings
import astropy
import hashlib

warnings.filterwarnings("ignore", module='astropy.io.votable.tree')
warnings.filterwarnings("ignore", category=astropy.table.TableReplaceWarning)

def hash(s):
    return hashlib.md5(s.encode()).hexdigest()[-6:]

def ic(x):
    # return
    print(str(x))

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



def get_gaiadr2_stars_new(ra, dec, distance=1, mag_min=5, mag_max=17):
    c = coord.SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
# coord = SkyCoord(ra=280, dec=-60, unit=(u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(c, u.Quantity(1.0, u.arcmin))
    r = j.get_results().to_pandas()
    r = r[(r['phot_g_mean_mag'] < mag_min) & (r['phot_g_mean_mag'] < mag_max)]    
    if len(r) == 0:
        return []
    r['star_mag'] = r['phot_g_mean_mag']
    return r

def get_gaiadr2_stars(ra, dec, distance=1, mag_min=18, mag_max=20):
#     print("Gaia DR2", mag_min, mag_max, distance)
    c = coord.SkyCoord(ra * u.deg, dec * u.deg)
    result = esasky.ESASky.query_region_catalogs(c, distance * u.arcmin,
                                                 "Gaia DR2")
    if len(result) == 0:
        return []
    result = result[0].to_pandas()
    stars = result[(result['phot_g_mean_mag'] > mag_min)
                   & (result['phot_g_mean_mag'] < mag_max)]
    stars['star_mag'] = tars['phot_g_mean_mag']
    return stars


def get_gs_stars(ra, dec, distance=1, mag_min=18, mag_max=20):
    c = coord.SkyCoord(ra * u.deg, dec * u.deg)
    my_catname = 'Guide Star Catalog 2.3 Cone Search 1'
    result = conesearch.conesearch(
        c, distance * u.arcmin, catalog_db=my_catname)
#     result = result.to_table().to_pandas()
    result = result.to_pandas()
    # print(len(result), 'objects')
    stars = result[(result.Mag >= mag_min) & (result.Mag <= mag_max) &
                   (result['class'] == 0)].copy()
    stars['star_mag'] = stars['Mag']
    return stars



def do_stars(lens, ra, dec, distance=1, mag_min=10, mag_max=14):
    try:
        objid = lens['name']
    except:
        objid = str(lens.name)

    stars = get_gaiadr2_stars_new(
        ra, dec, distance=distance, mag_min=mag_min, mag_max=mag_max)

    if len(stars) == 0:
        stars = get_gs_stars(ra, dec, distance=distance, mag_min=mag_min, mag_max=mag_max)
        if len(stars) == 0:
            # Can't find anything
            ic("No stars found")
            return None

    # for key in ["Mag", "phot_g_mean_mag", "VMag"]:
        # if key in stars.columns:
            # star_mag = stars.loc[closest, key]
            # break


    return stars


def update_cat(cat, i, stars):
    if stars is None or len(stars) == 0:
        cat.loc[i, 'star_mag'] = 99
        cat.loc[i, 'star_ra'] = -1000
        cat.loc[i, 'star_dec'] = -1000
        cat.loc[i, 'star_distance'] = 999999
        cat.loc[i, 'stars'] = 0
        return

    closest, dist = get_closest(stars, cat.loc[i, 'ra'], cat.loc[i, 'dec'])
    cat.loc[i, 'star_mag'] = stars.loc[closest, 'star_mag']
    cat.loc[i, 'star_ra'] = stars.loc[closest, 'ra']
    cat.loc[i, 'star_dec'] = stars.loc[closest, 'dec']
    cat.loc[i, 'star_distance'] = stars.loc[closest, 'distance_to_lens']
    cat.loc[i, 'stars'] = len(stars)


def get_closest(stars, ra, dec):
    print(ra, dec)
    closest = None
    cdistance = 1e9
    stars['distance_to_lens'] = 1000
    for id_, row in stars.iterrows():
        d = distance_sky((ra, dec), (row['ra'], row['dec']))
        stars.loc[id_, 'distance_to_lens'] = d
        if d < cdistance:
            cdistance = d
            closest = id_
    return closest, cdistance


def count_stars(cat_to_check, mag=16, distance=1, mag_min=10, use_cache=True):
    x = 0
    results = []
    for i, row in cat_to_check.iterrows():
        try:
            name = row['name']
        except:
            name = str(i)

        cache_check = f"# {name} {mag} {distance}"
        cache_hash = hash(cache_check)

        if not os.path.exists("stars"):
            os.mkdir("stars")
        fname = f"stars/{name}_stars_{cache_hash}.csv"
        ic(f"{i} {name} {cache_check} - {fname}")
        if use_cache and os.path.exists(fname):
            ic("Cached:" )
            ic(cache_check)
            ic(ct.fal(fname, skip_comments=False)[0] + "<<")
            # ic(cache_check == ct.fal(fname)[0])
            if cache_check == ct.fal(fname, skip_comments=False)[0]:
                ic("YES")
                stars = pd.read_csv(fname, comment="#")
                update_cat(cat_to_check, i, stars)
                continue

        stars = do_stars(dict(name=row['name']), row['ra'], row['dec'], distance=distance, mag_max=mag, mag_min=mag_min)
        update_cat(cat_to_check, i, stars)
        if stars is None:
            stars = pd.DataFrame()
        with open(fname, "w") as f:
            f.write(f"{cache_check}\n")
            stars.to_csv(f) 
    return
