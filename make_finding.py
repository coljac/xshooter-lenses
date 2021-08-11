import sys
import os
import time
import pandas as pd
import numpy as np
import coltools as ct
import pyds9
import pyregion
from PIL import Image
from PIL import Image, ImageDraw, ImageFont
from pyregion import Shape
import glob
import argparse

# fitsdir =  "/home/coljac/dev/science/slits/finding/fits"
# regions_dir =  "/home/coljac/dev/science/slits/regions/fixed3"
# regions_dir =  "/home/coljac/work/proposals/xshooter/2020A/targets/xshooter/region_files"
runid = "105.20KF.001"
piname = "S. Lopez"
fnt = ImageFont.truetype('/home/coljac/.local/share/fonts/arial.ttf', 24)
# outputdir = "/home/coljac/work/proposals/xshooter/2020A/targets/xshooter/finding_charts/"
# outputdir = "/home/coljac/dev/science/slits/finding/outputs"

def get_star_location(regions, index=1):
    x = 0
    for shape in regions:
        if shape.name == "line":
            x += 1
            if x == index:
                return shape.coord_list[0:2]
    return None

def get_slit_location(regions, index=1):
    x = 0
    for shape in regions:
        if shape.name == "box":
            x += 1
            if x == index:
                return shape.coord_list
    return None

def make_chart(name, output_file=None, outputdir=".", target_type="lens", scale="mode 99", size=768, fitsdir="fits", 
        regions_dir="."):

    # starmags = pd.read_csv("/home/coljac/Dropbox/Data/Development/science/slits/starmags.csv")
    # this_mag = starmags[starmags['name'] == name]
    # if len(this_mag) != 1:
        # print(len(this_mag), "error", name)
        # return
    # starmag = starmags[starmags['name'] == name].iloc[0]['starmag']

    d = pyds9.DS9()  # will open a new ds9 window or connect to an existing one
    d.set("frame delete all")
    d.set("frame new")
    d.set(f"width {size}")
    d.set(f"height {size}")
    
    fits = f"{fitsdir}/1024/{name}_g.fits"
    region_file = f"{regions_dir}/{name}.reg"

    if not os.path.exists(region_file):
        print(f"Region file for {name} is not ready. Skipping...")
        return
    regions = pyregion.open(region_file)
    slit_locations = {}
    
    star_location = get_star_location(regions)
    offset_coords = ct.decimal_to_hours(*star_location)
    offset_coords = f"{offset_coords[0].replace(' ', ':')} {offset_coords[1].replace(' ', ':')}"
    slit_locations['lens'] = get_slit_location(regions)
    slit_locations['source'] = get_slit_location(regions, index=3)
    
    d.set(f"file {fits}")
    
    d.set(f"scale {scale}")
    d.set("cmap invert yes")
    d.set("zoom to 1.3")
    d.set(f"pan to {star_location[0]} {star_location[1]} wcs")
    
    # Slit
    slit_location = slit_locations[target_type]

    linestr = f"fk5; line({star_location[0]} {star_location[1]} {slit_location[0]} {slit_location[1]}) # color=red width=2"
    slitstr = f"fk5; box({slit_location[0]} {slit_location[1]} 1.200\" 11.000\" {slit_location[4]} # color=red width=2"
    d.set("regions", linestr)
    d.set("regions", slitstr)

    # Decorations
    curr_x, curr_y = [float(x) for x in d.get("pan").split(" ")]
    compass_str = f"# compass({curr_x+150},{curr_y+150},10.520\") compass=fk5 {{N}} {{E}} 1 1 color=magenta   width=2 font=\"helvetica 24 normal roman\" fixed=1"
    curr_y += 100
    curr_x -= 70
    scale_str = f"line({curr_x},{curr_y},{curr_x+114},{curr_y}) # line=1 1 color=magenta width=3 font=\"helvetica 20 normal roman\" text={{30 arcsec}}"
    d.set("regions", compass_str)
    d.set("regions", scale_str)

    output_file = f"{outputdir}/finding_{name}_{target_type}.jpg"
    d.set(f"saveimage jpeg {output_file} 100")
    img = Image.open(output_file)
    w, h = img.size
    draw = ImageDraw.Draw(img)
    draw.rectangle([(0, h-47), (w, h)], fill="white")
    draw.text((0, h-40), f"runid: {runid}  PI: {piname}  Target: {name} ({target_type})  ", font=fnt, fill="black")
    # draw.text((0, h-160), f" starmag: {starmag: .1f}", font=fnt, fill="magenta")
    draw.text((0, h-200), f" λ=4000-5600", font=fnt, fill="magenta")
    draw.text((0, 50), f"Offset star coord: {offset_coords}", font=fnt, fill="magenta")
    d.set("width 256")
    d.set("height 256")
    d.set("zoom to 4")
    d.set(f"pan to {slit_location[0]} {slit_location[1]} wcs")
    d.set(f"saveimage jpeg /tmp/zoom.jpg 100")
    imz = Image.open("/tmp/zoom.jpg")
    imz = imz.crop((0, 0, imz.size[0], imz.size[1] - 47))
    img.paste(imz, (w-256, h-256-47))
    draw.line((w-256, h-256-47, w, h-256-47), fill="yellow", width=3)
    draw.line((w-256, h-256-47, w-256, h-47), fill="yellow", width=3)
    draw.line((w-256, h-47, w, h-47), fill="yellow", width=3)
    draw.line((w-3, h-47, w-3, h-256-47), fill="yellow", width=3)
    
    img.save(output_file)

def omain(argv):
    if len(argv) == 0:
        print("make_finding.py [-x <scale mode>] [-s <img size>] <ids as string or filename>")
        return

    if pyds9.ds9_targets() is None:
        print("No ds9, start it yourself!")
        return

    scale = "mode 99"
    scale = ct.argv(argv, "-x", "mode 99")
    size = int(ct.argv(argv, "-s", 768))

    ids = argv.pop(0)
    if "," in ids:
        ids = ids.split(",")
    else:
        ids = ct.fal(ids)

    for i in ids:
        if len(i) == 0:
            break
        for target in ['lens', 'source']:
            make_chart(i, target=target, 
                   scale=scale, 
                   size = size,
                   outputdir=outputdir)

def main(argv):
    parser = argparse.ArgumentParser(description="XShooter slits for lenses")
    parser.add_argument("--scale-mode", "-x", type=str, default="mode 99", help="Scale mode")
    parser.add_argument("--image-size", "-s", type=int, default=1024, help="Image size")
    parser.add_argument("--targets", "-t", type=str, required=False, help="Target names (list or file)")
    parser.add_argument("--output-dir", "-o", type=str, default=".", help="Target names (list or file)")
    parser.add_argument("--fits-dir", "-f", type=str, default="fits", help="Location of FITS files")
    args = parser.parse_args()

    if args.targets is None:
        targets = [x.replace(".reg", "") for x in glob.glob("*.reg")]
    else:
        if os.path.exists(args.targets):
            targets = ct.fal(args.targets)
        else:
            targets = args.targets.split(",")
    
    for target in targets:
        for target_type in ['lens', 'source']:
            print(f"Making chart for {target}, {target_type}")
            make_chart(target,
                    target_type=target_type,
                    scale=args.scale_mode,
                    size=args.image_size,
                    fitsdir=args.fits_dir,
                    outputdir=args.output_dir)

if __name__ == "__main__":
    main(sys.argv[1:])
