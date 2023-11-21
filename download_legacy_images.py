#Code for creating images from Legacy Survey DR9
#Written by Amanda Lopes

import argparse
import wget
from pathlib import Path
from astropy.table import Table
import os

parser = argparse.ArgumentParser(
    description="""Download images from Legacy""")
parser.add_argument("table", type=str,
                    default="table with the sources",
                    help="Name of table, taken the prefix ")

parser.add_argument("--Object", type=str,
                    default=None,
                    help="Id object of a given source ")
parser.add_argument("--legacy", action="store_true",help="make legacy images")
parser.add_argument("--DR", type=str,
                    default='dr9',
                    help="Data release to be used to obtain the images. Default: 'dr9'")
parser.add_argument("--bands", type=str,
                    default='grz',
                    help="Filters to be used in the combined image. Default: 'grz'")
parser.add_argument("--pixscale",type=str,
                    default='0.262',
                    help="Pixel scale of the images. Default: 0.262")

args = parser.parse_args()
file_ = args.table + ".ecsv"

try:
    data = Table.read(file_, format="ascii.ecsv")
except FileNotFoundError:
    file_ = args.table + ".csv"
    data = Table.read(file_, format="ascii")

if args.Object is not None:
    Object_ = str(args.Object)
    mask = np.array([source in Object_ for source in data["ID"]])
    data = data[mask]
else:
    data = data

def download_legacy(data,outpath):
    columns = ["RA","DEC","ID","radius"]
    for tab in data:
        try:
            ra = tab["RA"]
            dec = tab["DEC"]
            radii = tab["radius"]
            Name = tab["ID"]
        except KeyError as ke:
            print('No column name {} was found. Please check your table.'.format(ke))
            os._exit(0)

        print(Name)
        url = f"https://www.legacysurvey.org/viewer/jpeg-cutout?ra={ra}&dec={dec}&size={radii}&layer=ls-{args.DR}&pixscale={args.pixscale}&bands={args.bands}"
        wget.download(url, '{}/{}_{}pix.jpeg'.format(outpath,Name, radii))

if args.legacy==True:
    dir_output = Path(".")
    path_legacy = os.path.join(dir_output,'legacy_color_images')
    if os.path.exists(path_legacy)==False:
        os.mkdir(path_legacy)
        print("Directory '{}' created!".format(path_legacy))
    else:
        print("Directory '{}' already exists!".format(path_legacy))
    download_legacy(data,path_legacy)
    print()
    print("Download from Legacy Survey finished!")