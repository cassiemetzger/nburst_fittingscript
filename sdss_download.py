import numpy as np
import astropy 
from astropy.io import fits, ascii 
from astropy.table import Table, Column
import subprocess
import os
import glob
import argparse

SDSS_DESTINATION = '/Users/f007znp/Research/processed/SDSS_BOSS/'
def download_spectra(table): 
    for source in table: 
        if(source['SDSS_DR'] == '19.0'): 
            url = f"https://data.sdss.org/sas/dr19/spectro/sdss/redux/{source['sdss run']}/spectra/full/{source['sdss field']:06d}/{source['sdss mjd']}/{source['source file']}"
        else: 
            url = f"https://data.sdss.org/sas/dr17/sdss/spectro/redux/{source['sdss run']}/spectra/full/{source['sdss plate']:04d}/{source['source file']}"
        os.system(f"wget -P {SDSS_DESTINATION} {url}")
    

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str)
    args = parser.parse_args()
    hdu = fits.open(args.input)
    input = hdu[1].data
    m = input['survey'] == 'sdss'
    sdss_sources = input[m]
    sdss_files = glob.glob('/Users/f007znp/Research/processed/SDSS_BOSS/*')
    search_needed = [True if f"/Users/f007znp/Research/processed/SDSS_BOSS/{file}" not in sdss_files else False for file in sdss_sources['source file']]
    target_sources = sdss_sources[search_needed]
    download_spectra(target_sources)







if __name__ == "__main__":
    main()
