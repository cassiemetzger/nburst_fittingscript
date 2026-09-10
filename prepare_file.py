import numpy as np
import matplotlib.pyplot as plt
import astropy 
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.table import Table, Column
from astropy.cosmology import Planck15 as cosmo
import os
import requests
import numpy.ma as ma
import glob
import argparse
import pandas as pd 
from tqdm import tqdm

hdu = fits.open("./data/zall-tilecumulative-iron.fits")
desi=hdu[1].data
hdu=fits.open('./data/spAll-v6_1_3.fits')
dr19_list = hdu[1].data
dr19_specobj = [s.strip() for s in dr19_list['SPECOBJID'] ]
hdu=fits.open('./data/specObj-dr17.fits')
dr17_list = hdu[1].data
dr17_specobj = [s.strip() for s in dr17_list['SPECOBJID'] ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="The file you want to process")
    parser.add_argument("output", type=str, help="The desired output path")
    args = parser.parse_args()
    hdu = fits.open(args.input) 
    input_tbl = Table(hdu[1].data)
    desi_lookup = {tid: i for i, tid in enumerate(desi['TARGETID'])}
    dr17_lookup = {sid: i for i, sid in enumerate(np.asarray(dr17_specobj))}
    dr19_lookup = {sid: i for i, sid in enumerate(np.asarray(dr19_specobj))}

    survey = []
    runlist = []
    sourcefile = []
    tilelist = []
    lastnightlist = []
    fiberlist = []
    platelist = []
    fieldlist = []
    mjdlist = []
    fiberidlist = []
    fiberflux_z = []

    for i, source in enumerate(input_tbl):
        if source['DESI_TARGETID'] != '':
            survey.append('desi')
            idx = desi_lookup[int(source['DESI_TARGETID'])]
            tile = desi['TILEID'][idx]
            lastnight = desi['LASTNIGHT'][idx]
            petal = desi['PETAL_LOC'][idx]
            fiber = desi['FIBER'][idx]
            tilelist.append(tile)
            lastnightlist.append(lastnight)
            fiberlist.append(fiber)
            sourcefile.append(f"desiSp1d_{tile}-thru{lastnight}-{fiber}.fits")
            fiberflux_z.append(desi['FIBERFLUX_Z'][idx])
            runlist.append('none')
            platelist.append(9999)
            mjdlist.append(9999)
            fiberidlist.append(9999)
            fieldlist.append(9999)
            

        elif source['SDSS_SPECOBJID'] != '':
            survey.append('sdss')
            tilelist.append(9999)
            lastnightlist.append(9999)
            fiberlist.append(9999)
            fiberflux_z.append(9999)

            if source['SDSS_DR'] == '17.0':
                idx = dr17_lookup[source['SDSS_SPECOBJID']]
                plate = dr17_list['PLATE'][idx]
                runlist.append(dr17_list['RUN2D'][idx])
                fieldlist.append(9999)
                mjdlist.append(dr17_list['MJD'][idx])
                fiberidlist.append(dr17_list['FIBERID'][idx])
                sourcefile.append(f"spec-{plate:04d}-{dr17_list['MJD'][idx]}-{dr17_list['FIBERID'][idx]:04d}.fits")
                platelist.append(plate)

            elif source['SDSS_DR'] == '19.0':
                idx = dr19_lookup[source['SDSS_SPECOBJID']]
                runlist.append(dr19_list['RUN2D'][idx])
                sourcefile.append(dr19_list['SPEC_FILE'][idx])
                platelist.append(9999)
                fieldlist.append(dr19_list['FIELD'][idx])
                mjdlist.append(dr19_list['MJD'][idx])
                fiberidlist.append(dr19_list['FIBERID_LIST'][idx][0])

        else:
            survey.append('none')
            runlist.append('none')
            platelist.append(9999)
            fieldlist.append(9999)
            mjdlist.append(9999)
            fiberidlist.append(9999)
            tilelist.append(9999)
            lastnightlist.append(9999)
            fiberlist.append(9999)
            sourcefile.append('none')
            fiberflux_z.append(9999)
    input_tbl['survey'] = survey
    input_tbl['sdss run'] = runlist 
    input_tbl['sdss plate'] = platelist
    input_tbl['sdss field'] = fieldlist
    input_tbl['sdss mjd'] = mjdlist 
    input_tbl['sdss fiberid'] = fiberidlist
    input_tbl['desi tile'] = tilelist
    input_tbl['desi lastnight'] = lastnightlist
    input_tbl['desi fiberid'] = fiberlist
    input_tbl['desi fiberflux_z'] = fiberflux_z
    input_tbl['source file'] = sourcefile
    input_tbl.write(f'{args.output}', format = 'fits', overwrite = True)

if __name__ == "__main__":
    main()
    

