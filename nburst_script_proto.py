import numpy as np
import astropy 
from astropy.io import fits, ascii 
from astropy.table import Table, Column, vstack
import subprocess
import os
import glob
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from matplotlib.gridspec import GridSpec
import subprocess
from sys import stderr
import tempfile
import os
import numpy as np
import astropy 
from astropy.io import fits, ascii 
from astropy.table import Table, Column
import subprocess
import os
import glob
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import re
from tqdm import tqdm
import argparse

sdss_filepath = '/Users/f007znp/Research/processed/SDSS_BOSS'
desi_filepath = '/Users/f007znp/Research/processed/DESI'
output_location = '/Users/f007znp/Research/IMBH/nburst_fittingscript/outputs'
fitting_output_location_sdss = '/Users/f007znp/Research/pro/trial4/SDSS_BOSS'
fitting_output_location_desi = '/Users/f007znp/Research/pro/trial4/DESI'

class File:
    def __init__(self, filelist: list, outputpath: str):
        self.filelist = filelist
        self.outputpath = outputpath

    def write_output(self):
        with open(self.outputpath, 'w') as f:
            f.writelines(f"{file}\n" for file in self.filelist)

sp = ascii.read('./nburst_fittingscript/sdss_softening_param.txt')

def sdss_flux_mag(f, band): 
    idx = np.where(sp['filter']==band)[0][0]
    b = sp['b'][idx]
    return -2.5/np.log10(10) * (np.asinh((f/1E9)/(2*b)) + np.log10(b))
def desi_flux_mag(f): 
    return 22.5-2.5*np.log10(f)

def run_gdl_script(gdl_code, workdir=None, timeout=600, total = None, pbar_desc = "Fitting spectra"):
    full_code = "start.pro\n .r /Users/f007znp/Research/nbursts/idl/ppxf_nbursts_c.pro\n" + gdl_code + "\nexit\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pro', delete=False) as f:
        f.write(full_code)
        script_path = f.name
    process_re = re.compile(r"Processing\s+(\d+)\s*/\s*(\d+)")
    try:
        process = subprocess.Popen(
            ["gdl", script_path],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, cwd=workdir, bufsize=1
        )
        lines = []
        pbar = tqdm(total=total, desc=pbar_desc, unit="obj")
        last_i = 0
        for line in process.stdout:
            print(line, end="")   
            lines.append(line)
            m = process_re.search(line)
            if m:
                i, imax = int(m.group(1)), int(m.group(2))
                if pbar.total is None:
                    pbar.total = imax + 1
                pbar.update(i - last_i)
                last_i = i
        process.wait(timeout=timeout)
        pbar.close()
        return "".join(lines), ""
    finally:
        os.remove(script_path)

def sort_redshift(source, zcutoff):
    m = source['Rec_Z']> zcutoff 
    highz = source[m]
    lowz = source[~m]
    return lowz, highz 


def sdss_lrg_definition(source): 
    sdss_i = sdss_flux_mag(source['SDSS_CALIBFLUX_i'], 'i')
    sdss_z = sdss_flux_mag(source['SDSS_CALIBFLUX_z'], 'z')
    sdss_r = sdss_flux_mag(source['SDSS_CALIBFLUX_r'], 'r')
    lrg_izw = (sdss_i - sdss_z > 0.7) & (sdss_i - source['WISE_w1mpro'] > (2.143)*(sdss_i -sdss_z) - 0.2) & (sdss_z < 19.95) & (sdss_i > 19.9)
    lrg_riw = (sdss_r - sdss_i > 0.98) & (sdss_r - source['WISE_w1mpro'] > 2*(sdss_r - sdss_i)) & (sdss_i - sdss_z > 0.625) & (sdss_z < 19.95) & (sdss_i > 19.9)
    lrg = (lrg_izw) | (lrg_riw)
    return lrg 
def desi_lrg_definition(source): 
    desi_g = desi_flux_mag(source['DESI_FLUX_G'])
    desi_r = desi_flux_mag(source['DESI_FLUX_R'])
    desi_z = desi_flux_mag(source['DESI_FLUX_Z'])
    onea = (desi_z - source['WISE_w1mpro'] > 0.8 * (desi_r - desi_z) - 0.6)
    oneb = ((desi_g - source['WISE_w1mpro'] > 2.6) & (desi_g - desi_r > 1.4)) | (desi_r - source['WISE_w1mpro'] > 1.8)
    onec = (desi_r - desi_z > (desi_z - 16.83) * 0.45) & (desi_r - desi_z > (desi_z - 13.80) * 0.19) 
    oned = desi_r - desi_z > 0.7 
    onee = desi_flux_mag(source['desi fiberflux_z']) < 21.5 
    lrg = onea & oneb & onec & oned & onee 
    return lrg 

def compute_bic(hdu, k, target, window = 30, mask=None):
    idx = np.where(hdu[2].data["LINE_ID"][0] == target)[0][0]
    wavelength = hdu[2].data['WAVE'][0][idx]
    spec = hdu[1].data 
    wave = spec['WAVE'][0]
    flux = spec['FLUX'][0]
    fit_total = spec['FIT'][0]
    error = spec['ERROR'][0]
    mask = (wave > wavelength - window) & (wave < wavelength + window)
    if mask is not None:
        flux, fit_total, error = flux[mask], fit_total[mask], error[mask]
    chi2 = np.sum(((flux - fit_total) / error) ** 2)
    n = len(flux)
    bic = chi2 + k * np.log(n)
    return bic, chi2

class Fitting: 
    def __init__(self, file:str, survey:str, fit_stellarpop: bool, nl_sig_limit:float):
        self.file = file 
        self.survey = survey
        self.fit_stellarpop = fit_stellarpop 
        self.nl_sig_limit = nl_sig_limit 
        if self.survey == 'sdss': 
            self.outpath = fitting_output_location_sdss
        else: 
            self.outpath = fitting_output_location_desi
    def fit_nl(self): 
        if(self.fit_stellarpop): 
            gdl_code = (
                    f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=2,emexcl=0,"
                    f"emlt1=[1],start=[0,100,0,80,3000,-1.2],/force_sigma,"
                    f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4, siglimits=[200, {self.nl_sig_limit}],"
                    "path_ssp='/Users/f007znp/Research/stellar_templates/XSL/Kroupa/',"
                    f"prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/nl_'"
                    )
        else: 
             gdl_code = (
                        f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=2,emexcl=0,"
                        f"emlt1=[1],start=[0,100,0,80,3000,-1.2],/force_sigma,/disable_stpop,"
                        f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4,siglimits=[200, {self.nl_sig_limit}],"
                        f"outpath='{self.outpath}/nl_'"
                        )
        run_gdl_script(gdl_code, workdir="/Users/f007znp/Research/pro",total=1, pbar_desc = "Fitting NL")
    def fit_bl(self, sig_max): 
        if(self.fit_stellarpop): 
            gdl_code = (
                f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=3,emexcl=0,"
                f"emlt1=[1],emlt2 = [2], start=[0,100,0,80,0, 400, 3000,-1.2],/force_sigma,"
                f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4,siglimits=[200, {self.nl_sig_limit},{sig_max}],"
                "path_ssp='/Users/f007znp/Research/stellar_templates/XSL/Kroupa/',"
                f"prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/bl_'"
                )
        else: 
            gdl_code = (
            f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=3,"
            f"emlt1=[1], emlt2=[2], start=[0,100,0,80,0,400, 3000,-1.2],/force_sigma,/disable_stpop,"
            f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4,siglimits=[200, {self.nl_sig_limit},{sig_max}],"
            f"outpath='{self.outpath}/bl_'"
            )
        run_gdl_script(gdl_code, workdir="/Users/f007znp/Research/pro",total=1, pbar_desc = "Fitting BL")
    def fit_outflow(self, classification, start_velo, sig_max_o, sig_max): 
        if(self.fit_stellarpop): 
            if(classification in (0,-1)): 
                gdl_code = (
                        f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=3,emexcl=0,"
                            f"emlt1=[1,2],start=[0,100,0,80,0,150,3000,-1.2],/force_sigma,"
                            f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4,siglimits=[200,150,{sig_max_o}],symlosvdid=[2],noh4h6losvdid=[2],addstart=[0,0,0,0,-{start_velo},0,0,0],"
                            f"path_ssp='/Users/f007znp/Research/stellar_templates/XSL/Kroupa/',"
                            f"prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/out_'"
                            )
            else: 
                gdl_code = ( f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=4,emexcl=0,"
                            f"emlt1=[1,2],emlt2 =[3],start=[0,100,0,80,0,150,0,400, 3000,-1.2],/force_sigma,addstart=[0,0,0,0,-{start_velo},0,0,0,0,0],"
                            f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4, siglimits=[200,150,{sig_max_o}, {sig_max}],symlosvdid=[2],noh4h6losvdid=[2],"
                            f"path_ssp='/Users/f007znp/Research/stellar_templates/XSL/Kroupa/',"
                            f"prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/out_'")
        else: 
            if(classification in (0,-1)): 
                gdl_code = ( f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=3,emexcl=0,"
                            f"emlt1=[1,2],start=[0,100,0,80,0,150,3000,-1.2],/force_sigma,/disable_stpop,addstart=[0,0,0,0,-{start_velo},0,0,0],"
                            f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4, siglimits=[200,150,{sig_max_o}],symlosvdid=[2],noh4h6losvdid=[2],"
                            f"outpath='{self.outpath}/out_'")
            else: 
                gdl_code = (f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=4,emexcl=0,"
                            f"emlt1=[1,2],emlt2 =[3],start=[0,100,0,80,0,150,0,400, 3000,-1.2],symlosvdid=[2],noh4h6losvdid=[2],"
                            f"/force_sigma,/disable_stpop,addstart=[0,0,0,0,-{start_velo},0,0,0,0,0],"
                            f"lammin=3700,lammax=9000,degree=-1,mdegree=5,moments=4, siglimits=[200,150,{sig_max_o},{sig_max}],"
                            "outpath='/Users/f007znp/Research/pro/SDSS_BOSS/scripttest/out_'")
        run_gdl_script(gdl_code, workdir="/Users/f007znp/Research/pro",total=1, pbar_desc = "Fitting outflow")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="The file you want to process")
    parser.add_argument("output", type=str, help="The desired output path")
    args = parser.parse_args()

    hdu = fits.open(args.input)
    input = Table(hdu[1].data)
    if not os.path.exists(args.output) or os.path.getsize(args.output) == 0:
        with open(args.output, 'w') as f:
            f.write("file, ra, dec, redshift,survey, dropped, continuum, classification, delta_BIC, outflow, fit_fail, hit_limit, nburst_file" + '\n')
    nl_sig_limit = 150 
    finished = False
    for source in input[:2]: 
        while finished == False: 
            file = source['source file'] 
            ra = source['ra'] 
            dec = source['dec'] 
            redshift = source['Rec_Z'] 
            survey = source['survey']
            continuum = 0 
            classification = 0 
            delta_BIC = 9999
            outflow = 0 
            fit_failed = 0 
            hit_limit = 0 
            nburst_file = 'none'
            if((((source['SDSS_SPECOBJID'] != '')and (source['SDSS_CLASS'] != 'STAR')) or ((source['DESI_TARGETID'] != '')and (source['DESI_SPECTYPE'] != 'STAR'))) and (source['Rec_Z'] > 0)): 
                dropped_source = 0 
                if((source['Rec_Z']>0.7) and (source['survey'] == 'desi') and (desi_lrg_definition(source) == False)): 
                    continuum = 0 
                else: 
                    continuum = 1
            else: 
                dropped_source = 1   
            if(dropped_source == 0): 
                if(survey == 'sdss'):
                    with open(f'{output_location}/source_file.txt', 'w') as f: 
                        f.write(f"{sdss_filepath}/{file}")
                else: 
                    with open(f'{output_location}/source_file.txt', 'w') as f: 
                        f.write(f"{desi_filepath}/{file}")
                if(continuum == 0): 
                    fit_stellarpop = False 
                else: 
                    fit_stellarpop = True 

                fitter = Fitting(f"{output_location}/source_file.txt", survey, fit_stellarpop, nl_sig_limit)
                fitter.fit_nl()
                f = ascii.read(f"{output_location}/source_file.txt", format = 'no_header')[0][0]
                f = f.split("/")[6].split("-")
                one = f[1]
                two = f[2]
                three = f[3]
                if(source['survey'] == 'sdss'): 
                    nl_nburst_file = glob.glob(f"{fitting_output_location_sdss}/nl_/nbursts_sdss_*{one.strip("0")}*_*{two.strip("0")}*_*{three.strip("0")}*")
                else: 
                    nl_nburst_file = glob.glob(f"{fitting_output_location_desi}/nl_/nbursts_desi_*{one.strip("0")}*_*{two.strip("0")}*_*{three.strip("0")}*")
                hdu_nl = fits.open(nl_nburst_file[0])
                nl_sig = hdu_nl[1].data['SIG'][0][-1]
                bl_sig_min = 1.5*nl_sig
                sig_max = 5*nl_sig
                fitter.fit_bl(sig_max)
                if(source['survey'] == 'sdss'): 
                    bl_nburst_file = glob.glob(f"{fitting_output_location_sdss}/bl_/nbursts_sdss_*{one.strip("0")}*_*{two.strip("0")}*_*{three.strip("0")}*")
                else: 
                    bl_nburst_file = glob.glob(f"{fitting_output_location_desi}/bl_/nbursts_desi_*{one.strip("0")}*_*{two.strip("0")}*_*{three.strip("0")}*")
                n = hdu_nl[0].header['NWLFIT']
                dof = hdu_nl[0].header['DOF']
                k_nl= n-dof 
                hdu_bl = fits.open(bl_nburst_file[0])
                sig = hdu_bl[1].data['SIG'][0][2]
                if(sig < bl_sig_min or sig > 5*nl_sig): 
                    classification = 0 
                else: 
                    n = hdu_bl[0].header['NWLFIT']
                    dof = hdu_bl[0].header['DOF']
                    k_bl= n-dof 
                    lammin = hdu_nl[0].header['LAMMIN']
                    lammax = hdu_nl[0].header['LAMMAX']
                    if(lammax > 6562.8 and lammin< 6562.8):
                        target = "H alpha"
                        idx = np.where(hdu_nl[2].data['LINE_ID'][0] == target)[0][0]
                        nl_halpha_fl = hdu_nl[2].data['FLUX'][0][idx]
                        nl_halpha_err = hdu_nl[2].data['FLUX_ERR'][0][idx]
                        idx = np.where(hdu_bl[2].data['LINE_ID'][0] == target)[0][0]
                        bl_halpha_fl =  hdu_bl[2].data['FLUX'][0][idx]
                        bl_halpha_err = hdu_bl[2].data['FLUX_ERR'][0][idx]
                        if(np.isnan(nl_halpha_fl) or np.isnan(bl_halpha_fl)): 
                            fit_failed = 1 
                            bic_nl = np.nan
                            bic_bl = np.nan
                        else: 
                            if(nl_halpha_fl/nl_halpha_err >= 3): 
                                bic_nl, chi2_nl = compute_bic(hdu_nl, k_nl, target)
                            else: 
                                bic_nl = np.nan
                            if(bl_halpha_fl/bl_halpha_err >= 3):
                                bic_bl, chi2_bl = compute_bic(hdu_bl, k_bl, target)
                            else: bic_bl = np.nan
                    else: 
                        target = "Mg II] 2803"
                        idx = np.where(hdu_nl[2].data['LINE_ID'][0] == target)[0][0]
                        nl_mg_fl = hdu_nl[2].data['FLUX'][0][idx]
                        nl_mg_err = hdu_nl[2].data['FLUX_ERR'][0][idx]
                        idx = np.where(hdu_bl[2].data['LINE_ID'][0] == target)[0][0]
                        bl_mg_fl = hdu_bl[2].data['FLUX'][0][idx]
                        bl_mg_err = hdu_bl[2].data['FLUX_ERR'][0][idx]
                        if(np.isnan(nl_mg_fl) == False and np.isnan(bl_mg_fl) == False): 
                            if(nl_mg_fl/nl_mg_err >= 3): 
                                bic_nl_mgII, chi2_nl_mgII = compute_bic(hdu_nl, k_nl, target)
                            else: 
                                bic_nl_mgII = np.nan
                            if(bl_mg_fl/bl_mg_err >=3): 
                                bic_bl_mgII, chi2_bl_mgII = compute_bic(hdu_bl, k_bl, target)
                            else: 
                                bic_bl_mgII = np.nan
                        else: 
                            flag_fit = 1 
                            bic_nl_mgII = np.nan
                            bic_bl_mgII = np.nan
                        target = "H beta"
                        idx = np.where(hdu_nl[2].data['LINE_ID'][0] == target)[0][0]
                        nl_hbeta_fl = hdu_nl[2].data['FLUX'][0][idx]
                        nl_hbeta_err = hdu_nl[2].data['FLUX_ERR'][0][idx]
                        idx = np.where(hdu_bl[2].data['LINE_ID'][0] == target)[0][0]
                        bl_hbeta_fl = hdu_bl[2].data['FLUX'][0][idx]
                        bl_hbeta_err = hdu_bl[2].data['FLUX_ERR'][0][idx]
                        if(np.isnan(nl_hbeta_fl) == False and np.isnan(bl_hbeta_fl) == False): 
                            if(nl_hbeta_fl/nl_hbeta_err >= 3): 
                                bic_nl_hbeta, chi2_nl_hbeta = compute_bic(hdu_nl, k_nl, target)
                            else: 
                                bic_nl_hbeta = np.nan
                            if(bl_hbeta_fl/bl_hbeta_err >=3): 
                                bic_bl_hbeta, chi2_bl_hbeta = compute_bic(hdu_bl, k_bl, target)
                            else: 
                                bic_bl_hbeta = np.nan
                        else: 
                            flag_fit = 1 
                            bic_nl_hbeta = np.nan
                            bic_bl_hbeta = np.nan
                        if((np.isnan(nl_mg_fl) or np.isnan(bl_mg_fl) ) and (np.isnan(nl_hbeta_fl) or np.isnan(bl_hbeta_fl))): 
                            fit_failed = 1 
                        elif((np.isnan(nl_mg_fl) or np.isnan(bl_mg_fl)) and (np.isnan(nl_hbeta_fl) == False and np.isnan(bl_hbeta_fl) == False)):
                            bic_nl = nl_hbeta_fl 
                            bic_bl = bl_hbeta_fl
                        elif((np.isnan(nl_mg_fl) ==False and np.isnan(bl_mg_fl) == False) and (np.isnan(nl_hbeta_fl) or np.isnan(bl_hbeta_fl))):
                            bic_nl = nl_mg_fl 
                            bic_bl = bl_mg_fl
                        elif(np.isnan(nl_mg_fl) == False and np.isnan(bl_mg_fl) == False and np.isnan(nl_hbeta_fl) == False and np.isnan(bl_hbeta_fl) == False): 
                            bic_nl = np.mean([bic_nl_mgII, bic_nl_hbeta])
                            bic_bl = np.mean([bic_bl_mgII, bic_bl_hbeta])
                        else: 
                            bic_nl = np.nan
                            bic_bl = np.nan
                    if((np.isnan(bic_nl) == False) and (np.isnan(bic_bl)==False)): 
                        if(bic_nl < bic_bl): 
                            delta = bic_bl - bic_nl
                            if(delta > 6): 
                                classification = 0
                            else: 
                                classification = -1 
                        else: 
                            delta = bic_nl - bic_bl 
                            if(delta > 6): 
                                classification = 1
                            else: 
                                classification = -1
                    else: 
                        fit_failed = 1 
                if(fit_failed != 1): 
                    if(classification == 1): 
                        h3 = hdu_bl[1].data['H3'][0][1]
                        h3_err = hdu_bl[1].data['E_H3'][0][1]
                        sig = hdu_bl[1].data['SIG'][0][1]
                    else: 
                        h3 = hdu_nl[1].data['H3'][0][1]
                        h3_err = hdu_nl[1].data['E_H3'][0][1]
                        sig = hdu_bl[1].data['SIG'][0][1]
                    start_velo = 1.5*sig 
                    sig_max_o= 3*sig 
                    if(h3+h3_err < 0): 
                        fitter.fit_outflow(classification, start_velo, sig_max_o, sig_max)
                        if(source['survey'] == 'sdss'): 
                            out_nburst_file = glob.glob(f"{fitting_output_location_sdss}/out_/nbursts_sdss_*{one.strip("0")}*_*{two.strip("0")}*_*{three.strip("0")}*")
                        else: 
                            out_nburst_file = glob.glob(f"{fitting_output_location_desi}/out_/nbursts_desi_*{one.strip("0")}*_*{two.strip("0")}*_*{three.strip("0")}*")
                        hdu_out = fits.open(out_nburst_file[0])
                        n = hdu_out[0].header['NWLFIT']
                        dof = hdu_out[0].header['DOF']
                        k_out = n-dof
                        if(classification in (0,-1)): 
                            original = hdu_nl
                            k_original = k_nl
                        else: 
                            original = hdu_bl 
                            k_original = k_bl
                        lammin = hdu_out[0].header['LAMMIN']
                        lammax = hdu_out[0].header['LAMMAX']
                        if(lammin <= 5008 and lammax >= 5008): 
                            target = "[O III] 5008"
                            idx = np.where(original[2].data['LINE_ID'][0] ==target)[0][0]
                            fl = original[2].data['FLUX'][0][idx]
                            fl_err = original[2].data['FLUX_ERR'][0][idx]
                            if(fl/fl_err >= 3): 
                                bic_orig, chi2_orig = compute_bic(original, k_original, target)
                            else: 
                                bic_orig = np.nan 
                            idx = np.where(hdu_out[2].data['LINE_ID'][0] ==target)[0][0]
                            fl = hdu_out[2].data['FLUX'][0][idx]
                            fl_err = hdu_out[2].data['FLUX_ERR'][0][idx]
                            if(fl/fl_err >= 3): 
                                bic_out, chi2_out = compute_bic(hdu_out, k_out, target)
                            else: 
                                bic_out = np.nan
                            if(np.isnan(bic_orig) or np.isnan(bic_out)): 
                                fit_failed = 1 
                                outflow = 0
                            else: 
                                if(bic_orig > bic_out): 
                                    delta_out = bic_orig - bic_out 
                                else:
                                    delta_out = bic_out - bic_orig
                                if(delta_out > 6): 
                                    outflow = 1
                    if(classification in (0,-1) and outflow == 0): 
                        nburst_file = nl_nburst_file
                    if(classification == 1 and outflow == 0): 
                        nburst_file = bl_nburst_file
                    if(outflow ==1): 
                        nburst_file = out_nburst_file
                    if(nburst_file != 'none'):
                        print(nburst_file)
                        hdu = fits.open(nburst_file[0])
                        sig = hdu[1].data['SIG'][0][1][0]
                        if(np.isclose(sig, nl_sig_limit, atol=1e-1)): 
                            if(2*nl_sig_limit < 850):
                                nl_sig_limit= 2*nl_sig_limit
                            else: 
                                nl_sig_limit = 850
                                finished = True
                            hit_limit = 1 
                        else: 
                            finished = True 
                else: 
                    finished = True 
            else: 
                finished = True 

        line = f"{file}, {ra}, {dec}, {redshift}, {survey}, {continuum}, {classification}, {delta_BIC}, {outflow}, {fit_failed}, {hit_limit}, {nburst_file[0]}\n"
        with open(args.output, 'a') as f: 
            f.write(line)      



    """

    # remove sources without sdss or desi counterparts and stars 
    dropped_sources = np.asarray([0]*len(input))
    m = (((input['SDSS_SPECOBJID'] != '')& (input['SDSS_CLASS'] != 'STAR')) | ((input['DESI_TARGETID'] != '')& (input['DESI_SPECTYPE'] != 'STAR'))) & (input['Rec_Z'] > 0)
    dropped_sources[m] = 1 
    output_tbl['Dropped'] = dropped_sources

    survey_full = np.array(['none'] * len(input), dtype=object)
    survey_full[m & (input['DESI_TARGETID'] != '')] = 'desi'
    survey_full[m & (input['SDSS_SPECOBJID'] != '') & (survey_full != 'desi')] = 'sdss'
    output_tbl['survey'] = survey_full
    continuum = np.array([1] * len(input))
    redshift_sort_mask = (input['Rec_Z']>0.7) & (output_tbl['survey'] == 'desi') & (desi_lrg_definition(input) == False)
    continuum[redshift_sort_mask] = 0
    redshift_sort_mask = (input['Rec_Z']>0.7) & (output_tbl['survey'] == 'sdss') & (sdss_lrg_definition(input) == False)
    continuum[redshift_sort_mask] = 0
    output_tbl['continuum'] = continuum
    input = input[m]
    

    # sort by survey 
    m_desi = input['survey'] == 'desi'
    desilist = input[m_desi]
    m_sdss = input['survey'] == 'sdss'
    sdsslist = input[m_sdss]

    # sort by redshift for each survey 
    zcutoff = 0.7
    lowz_desi, highz_desi = sort_redshift(desilist, zcutoff)
    lowz_sdss, highz_sdss = sort_redshift(sdsslist, zcutoff)

    # check whether highz sources are lrgs 
    m = desi_lrg_definition(highz_desi)
    lowz_desi = vstack([lowz_desi, highz_desi[m]]) # add lrgs to lowz source list as they should be fit with stellar continuum 
    highz_desi = highz_desi[~m]
    m = sdss_lrg_definition(highz_sdss)
    lowz_sdss = vstack([lowz_sdss, highz_sdss[m]]) # add lrgs to lowz source list as they should be fit with stellar continuum 
    highz_sdss = highz_sdss[~m]
    
    
    # output into file lists to be fit 
    desi_filelist = [
            f"{desi_filepath}/{s['desi tile']}/{s['desi lastnight']}/1d/{s['source file']}"
            for s in lowz_desi]
    filelist = File(desi_filelist, f'{output_location}/desi_files_lowz.txt')
    filelist.write_output()
    desi_filelist = [
                f"{desi_filepath}/{s['desi tile']}/{s['desi lastnight']}/1d/{s['source file']}"
                for s in highz_desi]
    filelist = File(desi_filelist, f'{output_location}/desi_files_highz.txt')
    filelist.write_output()
    sdss_filelist = [
            f"{sdss_filepath}/{s['source file']}"
            for s in lowz_sdss]
    filelist = File(sdss_filelist, f'{output_location}/sdss_files_lowz.txt')
    filelist.write_output()
    sdss_filelist = [
                f"{sdss_filepath}/{s['source file']}"
                for s in highz_sdss]
    filelist = File(sdss_filelist, f'{output_location}/sdss_files_highz.txt')
    filelist.write_output()

    # fit 
    fitter = Fitting(f'{output_location}/sdss_files_lowz.txt', 'sdss', True, 150)
    fitter.run_full()


    ascii.write(output_tbl, args.output, format = 'csv', overwrite = True)





    m_desi = (input['survey'] == 'desi') & (input['DESI_SPECTYPE'] != 'STAR') & (input['Rec_Z'] > 0.7) 

    desi_filelist = [
        f"{desi_filepath}/{s['desi tile']}/{s['desi lastnight']}/1d/{s['source file']}"
        for s in input[m_desi]]
    filelist = File(desi_filelist, f'{output_location}/desi_files.txt')
    filelist.write_output()

    m_sdss = (input['survey'] == 'sdss') * (input['SDSS_CLASS'] != 'STAR') & (input['Rec_Z'] > 0)
    sdss_filelist = [
        f"{sdss_filepath}/{s['source file']}"
        for s in input[m_sdss]]
    filelist = File(sdss_filelist, f'{output_location}/sdss_files.txt')
    filelist.write_output()
    """

    




if __name__ == "__main__":
    main()