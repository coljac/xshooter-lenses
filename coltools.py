import sys
import os
import re
import numpy as np

# file-as-lines (list)
def fal(filename, skip_comments=True):
    with open(filename, "r") as f:
        lines = f.readlines()
    if skip_comments:
        return [l.strip() for l in lines if not l.startswith("#")]
    else:
        return [l.strip() for l in lines]


# File-as-string
def fas(filename):
    with open(filename, "r") as f:
        data = f.read()  # .replace('\n', '')
    return data


# file-as-dataframe
def fad(filename, format="csv", skip_comments=True):
    import pandas as pd
    if skip_comments:
        return pd.read_csv(filename, comment='#')
    else:
        return pd.read_csv(filename)


# save as dataframe
def sad(frame, filename):
    frame.to_csv(filename)


# File-as-table (astropy tables)
def fat(filename, format="csv"):
    from astropy.io import ascii
    return ascii.read(filename, format=format)


# Save-as-table (save table to file)
def sat(table, filename, format="csv"):
    from astropy.io import ascii
    ascii.write(table, filename, format=format)


# string to file
def strf(filename, string):
    with open(filename, "w") as f:
        f.write(string)


def printerr(string):
    sys.stderr.write(string + "\n")


def printc(string, color, nl=True, stream=sys.stdout):
    from termcolor import colored
    stream.write(colored(string, color))
    if nl:
        stream.write("\n")
    stream.flush()


def progbar(current, to, width=40, show=True, message=None, stderr=False):
    percent = float(current) / float(to)
    length = int(width * percent)
    if show:
        count = " (%d/%d)    " % (current, to)
    else:
        count = ""
    if message:
        count += message
    outstream = sys.stderr if stderr else sys.stdout
    outstream.write(("\r[" + ("#" * length) + " " * (width - length) +
                     "] %0d" % (percent * 100)) + "%" + count)
    outstream.flush()


def progfile(done, out_of=100):
    try:
        jobid = os.environ['PBS_JOBID']
        if jobid is not None and len(jobid) > 0:
            progfile = "/home/cjacobs/.jobs/progress.%s" % jobid.split(".")[0]
            strf(progfile, "%d/%d\n" % (done, out_of))
    except KeyError as e:
        pass


def list_files(dirname, extension, only_names=False):
    import glob
    # onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files = glob.glob(dirname + os.sep + "*.%s" % (extension))
    if only_names:
        return [f.split("/")[-1] for f in files]
    else:
        return files


def argb(argv, switch):
    if switch in argv:
        del argv[argv.index(switch)]
        return True
    return False


def argv(argv, switch, default=None):
    if switch in argv:
        val = argv[argv.index(switch) + 1]
        del argv[argv.index(switch) + 1]
        del argv[argv.index(switch)]
        return val
    return default


def write_with_lock(text, filename):
    import fcntl
    with open(filename, "a") as g:
        fcntl.flock(g, fcntl.LOCK_EX)
        g.write(text)
        fcntl.flock(g, fcntl.LOCK_UN)


def find_files(dirname, match=None):
    result = []
    for root, dirs, files in os.walk(dirname):
        for file_ in files:
            if match is None or match in file_:
                result.append(root + "/" + file_)
    return result


def decimal_to_hours(ra, dec, return_tuple=False, delim=" "):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    c = SkyCoord(ra, dec, unit="deg")
    if return_tuple:
        return [d for d in c.ra.hms] + [
            " " if c.dec.signed_dms[0] > 0 else " -"
        ] + [d for d in c.dec.signed_dms[1:]]
    fullstring = "%02d:%02d:%.4f %s%02d:%02d:%.4f" % (c.ra.hms + tuple(
        [" " if c.dec.signed_dms[0] > 0 else " -"]) + c.dec.signed_dms[1:])
    ra_hms = [d for d in c.ra.hms]
    dec_dms = [d for d in c.dec.signed_dms]
    fullstring = "%s %s %s_%s %s %s" % (str(int(ra_hms[0])).rjust(2, "0"), \
            str(int(ra_hms[1])).rjust(2, "0"), \
            '{:05.2f}'.format(ra_hms[2]), \
            ("-" if dec_dms[0] < 0 else "") + str(int(dec_dms[1])).rjust(2, "0"), \
            str(int(dec_dms[2])).rjust(2, "0"), '{:05.2f}'.format(dec_dms[3]))
    
    # (c.ra.hms + tuple([" " if c.dec.signed_dms[0] >0 else " -"]) + c.dec.signed_dms[1:])
    result = fullstring.split("_")
    result = [x.replace(" ", delim) for x in result]
    return result


def hours_to_decimal(ra, dec=""):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    c = SkyCoord("%s %s" % (ra, dec), unit=(u.hourangle, u.deg))
    return c.ra.value, c.dec.value

def distance_radec(ra, dec, ra2, dec2):
    coords1 = (ra, dec)
    coords2 = (ra2, dec2)
    offset_ra = 3600 * (coords2[0] - coords1[0]) * np.cos(coords2[1]*np.pi/180)
    offset_dec = 3600 * (coords2[1] - coords1[1])
    return np.sqrt(np.power(offset_ra, 2) + np.power(offset_dec, 2)) # ARCSEC

def get_sources_nearby(radec, catalog, distance=10, radial=True):  # arcsec
    """ catalog: pandas dataframe with ra, dec columns
        distance: Distance in arcsec"""

    ra, dec = radec
    distance = distance / 3600.
    if not radial:  # in a box
        ras = (ra - distance, ra + distance)
        decs = (dec - distance, dec + distance)
        sources = catalog[(catalog.ra > ras[0]) & (catalog.ra < ras[1]) &\
                (catalog.dec > decs[0]) & (catalog.dec < decs[1])]
    else:
        import numpy as np
        catras = catalog['ra']
        catdecs = catalog['dec']
        radist = catras - ra
        decdist = catdecs - dec
        distances2 = np.power(radist, 2) + np.power(decdist, 2)
        sources = catalog[distances2 <= distance**2]
    return sources


class Cobj(dict):
    def __getattr__(self, item):
        if item in self:
            return self[item]
        return None

    def __init__(self):
        pass

    def __dir__(self):
        return super().__dir__() + [str(k) for k in self.keys()]


def get_file_band(filename):
    band = None
    r = re.compile(r'([_-])([grizY])([_\.-])')
    try:
        m = re.search(r, filename)
        # m.group(0), m.group(1), m.group(2), m.group(3)
        band = m.group(2)
    except:
        pass
    return band


def get_matching_band(filename, band):
    r = re.compile(r'([_-])([grizY])([_\.-])')
    try:
        thisband = get_file_band(filename)
        if thisband is not None:
            a = re.sub(r, r"\1" + band + r"\3", filename)
            return a
    except:
        pass

    return None

def to_hours_name(ra, dec, prefix=""):
    d2h = decimal_to_hours(ra, dec)
    # hours_string = decimal_to_hours(ra, dec) #.split("  ")
    hours_string = d2h[0].split(" ")
    dec = d2h[1].split(" ")
    ra = hours_string
    ra = "".join([a.rjust(2,"0") for a in ra])[0:4]
    # dec = hours_string[1].split(" ")

    # ra = hours_string[0].split(":")
    # ra = "".join([a.rjust(2,"0") for a in ra])[0:4]
    # dec = hours_string[1].split(":")

    if not str(dec[0])[0] == "-":
        dec = ["+"] + dec
    else:
        dec = ["-"] + [dec[0][1:]] + dec[1:]
    dec = dec[0] + "".join([a.rjust(2,"0") for a in dec[1:3]])
    return prefix + ra+dec

def extract_pattern(string, regex):
    r = re.compile(regex)
    m = re.search(r, string)
    return m.group(1)

