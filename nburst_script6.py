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
            if(classification == 'nl'): 
                gdl_code = (
                        f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=3,emexcl=0,"
                            f"emlt1=[1,2],start=[0,100,0,80,0,150,3000,-1.2],/force_sigma,"
                            f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4,siglimits=[200,150,{sig_max_o}],symlosvdid=[2],noh4h6losvdid=[2],addstart=[0,0,0,0,-{start_velo},0,0,0],"
                            f"path_ssp='/Users/f007znp/Research/stellar_templates/XSL/Kroupa/',"
                            f"prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/out_'"
                            )
            elif(classification == 'bl'): 
                gdl_code = ( f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=4,emexcl=0,"
                            f"emlt1=[1,2],emlt2 =[3],start=[0,100,0,80,0,150,0,400, 3000,-1.2],/force_sigma,addstart=[0,0,0,0,-{start_velo},0,0,0,0,0],"
                            f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4, siglimits=[200,150,{sig_max_o}, {sig_max}],symlosvdid=[2],noh4h6losvdid=[2],"
                            f"path_ssp='/Users/f007znp/Research/stellar_templates/XSL/Kroupa/',"
                            f"prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/out_'")
        else: 
            if(classification == 'nl'): 
                gdl_code = ( f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=3,emexcl=0,"
                            f"emlt1=[1,2],start=[0,100,0,80,0,150,3000,-1.2],/force_sigma,/disable_stpop,addstart=[0,0,0,0,-{start_velo},0,0,0],"
                            f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4, siglimits=[200,150,{sig_max_o}],symlosvdid=[2],noh4h6losvdid=[2],"
                            f"outpath='{self.outpath}/out_'")
            elif(classification == 'bl'): 
                gdl_code = (f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=4,emexcl=0,"
                            f"emlt1=[1,2],emlt2 =[3],start=[0,100,0,80,0,150,0,400, 3000,-1.2],symlosvdid=[2],noh4h6losvdid=[2],"
                            f"/force_sigma,/disable_stpop,addstart=[0,0,0,0,-{start_velo},0,0,0,0,0],"
                            f"lammin=3700,lammax=9000,degree=-1,mdegree=5,moments=4, siglimits=[200,150,{sig_max_o},{sig_max}],"
                            "outpath='/Users/f007znp/Research/pro/SDSS_BOSS/scripttest/out_'")
        run_gdl_script(gdl_code, workdir="/Users/f007znp/Research/pro",total=1, pbar_desc = "Fitting outflow")
def divide_into_continuum_groups(table, definition): 
    m = table['Rec_Z'] > 0.7 
    highz = table[m]
    lowz = table[~m]
    if(definition == 'sdss'): 
        lrg = sdss_lrg_definition(highz)
    else: 
        lrg = desi_lrg_definition(highz)
    continuum = vstack([lowz, highz[lrg]])
    no_continuum = highz[~lrg]
    return continuum, no_continuum

def write_output_tbl(table): 
    tbl = Table(
        [table['name'], table['ra'], table['dec'], table['Rec_Z'], table['source file'], 
         table['continuum'], table['best fit'], table['outflow'], table['error message']],
        names=('name', 'ra', 'dec', 'Rec_Z', 'source_file', 'continuum', 'classification', 'outflow', 'error_message')
    )
    tbl.write(f"./nburst_fittingscript/outputs.txt", format='ascii', overwrite=True)

def retrieve_nburstfile(source, survey, spectype): 
    if(survey == 'desi'): 
            file_loc = fitting_output_location_desi.split("/")
            error_file = glob.glob(f"{file_loc[0]}/{file_loc[1]}/{file_loc[2]}/{file_loc[3]}/{file_loc[4]}/*{str(source['desi tile']).strip('0')}*_*{str(source['desi lastnight']).strip('0')}*_*{str(source['desi fiberid']).strip('0')}*") 
            if(len(error_file)>=1):
                file = 'None' 
                error_message = f"[ERROR] Fit for file {source['source file']} failed"
            else:
                file = glob.glob(f"{fitting_output_location_desi}/{spectype}/nbursts_desi_*{str(source['desi tile']).strip('0')}*_*{str(source['desi lastnight']).strip('0')}*_*{str(source['desi fiberid']).strip('0')}*") 
                error_message=''
                if(len(file) == 0): 
                    file = 'None'
                    error_message = f"[ERROR] Could not locate NBURSTS file associated with {source['source file']} "
    else: 
        if(source['dr'] == 17):
            file_loc = fitting_output_location_sdss.split("/")
            error_file = glob.glob(f"{file_loc[0]}/{file_loc[1]}/{file_loc[2]}/{file_loc[3]}/{file_loc[4]}/*{str(source['sdss plate']).strip('0')}*_*{str(source['sdss mjd']).strip('0')}*_*{str(source['sdss fiberid']).strip('0')}*") 
            if(len(error_file) >=1): 
                file = 'None'
                error_message = f"[ERROR] Fit for file {source['source file']} failed"
            else:
                file = glob.glob(f"{fitting_output_location_sdss}/{spectype}/nbursts_sdss_*{str(source['sdss plate']).strip('0')}*_*{str(source['sdss mjd']).strip('0')}*_*{str(source['sdss fiberid']).strip('0')}*")
                error_message=''
                if(len(file) == 0): 
                    file = 'None'
                    error_message = f"[ERROR] Could not locate NBURSTS file associated with {source['source file']} "
        else: 
            file_loc = fitting_output_location_sdss.split("/")
            sdss_file = source['source file'].split('-')
            error_file = glob.glob(f"{file_loc[0]}/{file_loc[1]}/{file_loc[2]}/{file_loc[3]}/{file_loc[4]}/*{sdss_file[1].strip('0')}*_*{sdss_file[2].strip('0')}*_*{sdss_file[3].split('.')[0].strip('0')}*")
            if(len(error_file) >=1): 
                file = 'None'
                error_message = f"[ERROR] Fit for file {source['source file']} failed"
            else:
                file = glob.glob(f"{fitting_output_location_sdss}/{spectype}/nbursts_sdss_*{sdss_file[1].strip('0')}*_*{sdss_file[2].strip('0')}*_*{sdss_file[3].split('.')[0].strip('0')}*")
                error_message = ''
                if(len(file) == 0): 
                    file = glob.glob(f"{fitting_output_location_sdss}/{spectype}/nbursts_sdss_*{sdss_file[1].strip('0')}*_*{sdss_file[2].strip('0')}*_*0000*")
                    error_message = ''
                    if(len(file) == 0): 
                        file = 'None'
                        error_message = f"[ERROR] Could not locate NBURSTS file associated with {source['source file']} "
    return file, error_message
def group_by_similarity(sig, table, tol=20.0):
    sig = np.asarray(sig, dtype=float)
    valid_mask = ~np.isnan(sig)
    valid_idx = np.where(valid_mask)[0]

    groups = []
    if len(valid_idx) > 0:
        order = valid_idx[np.argsort(sig[valid_idx])]
        current_group = [order[0]]
        for i in range(1, len(order)):
            idx_prev, idx_curr = order[i-1], order[i]
            if sig[idx_curr] - sig[idx_prev] <= tol:
                current_group.append(idx_curr)
            else:
                groups.append(current_group)
                current_group = [idx_curr]
        groups.append(current_group)

    nan_idx = np.where(~valid_mask)[0]
    if len(nan_idx) > 0:
        groups.append(list(nan_idx))

    return groups
def compare_fits(source, survey, target, file1path, file2path): 
    if(file1path != 'None' and file2path!= 'None'): 
        hdu_nl = fits.open(file1path[0])
        hdu_bl = fits.open(file2path[0])
        n = hdu_nl[0].header['NWLFIT']
        dof = hdu_nl[0].header['DOF']
        k_nl= n-dof
        n = hdu_bl[0].header['NWLFIT'] 
        dof = hdu_nl[0].header['DOF']
        k_bl = n-dof
        lammin = hdu_nl[0].header['LAMMIN']
        lammax = hdu_nl[0].header['LAMMAX']
        idx = np.where(hdu_nl[2].data['LINE_ID'][0] == target)[0][0]
        if(hdu_nl[2].data['WAVE'][0][idx] >= lammin and hdu_nl[2].data['WAVE'][0][idx] <= lammax): 
            idx = np.where(hdu_nl[2].data['LINE_ID'][0] == target)[0][0]
            nl_fl = hdu_nl[2].data['FLUX'][0][idx]
            nl_fl_err = hdu_nl[2].data['FLUX_ERR'][0][idx]
            idx = np.where(hdu_bl[2].data['LINE_ID'][0] == target)[0][0]
            bl_fl = hdu_bl[2].data['FLUX'][0][idx]
            bl_fl_err = hdu_bl[2].data['FLUX_ERR'][0][idx]
            if(np.isnan(nl_fl) or np.isnan(nl_fl_err)): 
                bic_nl = np.nan
            else: 
                if(nl_fl/nl_fl_err >= 3): 
                    bic_nl, chi2_nl = compute_bic(hdu_nl, k_nl, target)
                else: 
                    bic_nl = np.nan
            if(np.isnan(bl_fl) or np.isnan(bl_fl_err)): 
                bic_bl = np.nan
            else: 
                if(bl_fl/bl_fl_err >= 3): 
                    bic_bl, chi2_bl = compute_bic(hdu_bl, k_bl, target)
                else: 
                    bic_bl = np.nan
        else: 
            bic_nl = np.nan
            bic_bl = np.nan
    else: 
        bic_nl = np.nan
        bic_bl = np.nan
    return bic_nl, bic_bl 

def add_outflow(table, survey, continuum, classification, nl_sig_limit): 
    sig = []
    for source in table: 
        file, error_message = retrieve_nburstfile(source, survey, f"{classification}_")
        if file != 'None': 
            hdu = fits.open(file[0])
            sig.append(hdu[1].data['SIG'][0][1][0])
        else: 
            sig.append(np.nan)
    groups = group_by_similarity(sig, table)
    for i, idx_list in enumerate(groups):
        group = table[idx_list]
        med_sig = np.median(np.asarray(sig)[idx_list])
        bl_sig_max = 5 * med_sig 
        nl_offset = 1.5*med_sig 
        out_sig_max = 3*med_sig
        batch_path = f"{output_location}/{survey}_outflow_batch{i}.txt"
        for j, source in enumerate(group['source file']): 
            with open(batch_path, 'a') as f: 
                if survey == 'sdss': 
                    f.write(f"{sdss_filepath}/{source}\n")
                else: 
                    f.write(f"{desi_filepath}/{group['desi tile'][j]}/{group['desi lastnight'][j]}/1d/{source}\n")
        if os.path.exists(batch_path) and os.path.getsize(batch_path) > 0:
            batch_fit = Fitting(batch_path, survey, continuum, nl_sig_limit)
            batch_fit.fit_outflow(classification, nl_offset, out_sig_max,bl_sig_max)

def perform_fitting(table, survey, nl_sig_limit, bic_threshold): 
    if(f"{output_location}/{survey}_continuum.txt" in glob.glob(f"{output_location}/*")): 
        continuum_fit = Fitting(f"{output_location}/{survey}_continuum.txt", survey, True, nl_sig_limit)
        continuum_fit.fit_nl()
        m = table['continuum'] == True
        sig = []
        for source in table[m]: 
            file, error_message = retrieve_nburstfile(source, survey, 'nl_')
            if file != 'None': 
                hdu = fits.open(file[0])
                sig.append(hdu[1].data['SIG'][0][1][0])
            else: 
                sig.append(np.nan)
        groups = group_by_similarity(sig, table[m])
        subset = table[m]   # take this once, so you're not re-filtering every iteration
        for i, idx_list in enumerate(groups): 
            group = subset[idx_list]        # <-- index with ONE group's flat list, not all groups at once
            med_sig = np.median(np.asarray(sig)[idx_list])
            bl_sig_max = 5 * med_sig 
            batch_path = f"{output_location}/{survey}_bl_batch{i}.txt"
            for j, source in enumerate(group['source file']): 
                with open(batch_path, 'a') as f: 
                    if survey == 'sdss': 
                        f.write(f"{sdss_filepath}/{source}\n")
                    else: 
                        f.write(f"{desi_filepath}/{group['desi tile'][j]}/{group['desi lastnight'][j]}/1d/{source}\n")
            if os.path.exists(batch_path) and os.path.getsize(batch_path) > 0:
                batch_fit = Fitting(batch_path, survey, True, nl_sig_limit)
                batch_fit.fit_bl(bl_sig_max)

    if(f"{output_location}/{survey}_nocontinuum.txt" in glob.glob(f"{output_location}/*")): 
        nocontinuum_fit = Fitting(f"{output_location}/{survey}_nocontinuum.txt", survey, False, nl_sig_limit)
        nocontinuum_fit.fit_nl()
        sig = []
        for source in table[~m]: 
            file, error_message = retrieve_nburstfile(source, survey, 'nl_')
            if file != 'None': 
                hdu = fits.open(file[0])
                sig.append(hdu[1].data['SIG'][0][1][0])
            else: 
                sig.append(np.nan)
        groups = group_by_similarity(sig, table[~m])
        subset = table[m]   # take this once, so you're not re-filtering every iteration
        for i, idx_list in enumerate(groups): 
            group = subset[idx_list]        # <-- index with ONE group's flat list, not all groups at once
            med_sig = np.median(np.asarray(sig)[idx_list])
            bl_sig_max = 5 * med_sig 
            batch_path = f"{output_location}/{survey}_bl_batch{i}.txt"
            for j, source in enumerate(group['source file']): 
                with open(batch_path, 'a') as f: 
                    if survey == 'sdss': 
                        f.write(f"{sdss_filepath}/{source}\n")
                    else: 
                        f.write(f"{desi_filepath}/{group['desi tile'][j]}/{group['desi lastnight'][j]}/1d/{source}\n")
            if os.path.exists(batch_path) and os.path.getsize(batch_path) > 0:
                batch_fit = Fitting(batch_path, survey, True, nl_sig_limit)
                batch_fit.fit_bl(bl_sig_max)
    fit = []
    outflow_mask_nl_cont = [False]*len(table)
    outflow_mask_bl_cont = [False]*len(table)
    outflow_mask_nl_nocont = [False]*len(table)
    outflow_mask_bl_nocont = [False]*len(table)
    for i, source in enumerate(table): 
        file_nl, error_message = retrieve_nburstfile(source, source['survey'], 'nl_')
        file_bl, error_message = retrieve_nburstfile(source, source['survey'], 'bl_')
        if(file_nl != 'None' and file_bl!='None'): 
            if(source['Rec_Z'] < 0.371):
                bic_nl, bic_bl = compare_fits(source, survey, "H alpha", file_nl, file_bl)
                if(np.isnan(bic_nl) or np.isnan(bic_bl)): 
                    bic_nl, bic_bl = compare_fits(source, survey, "H beta", file_nl, file_bl)
                    if(np.isnan(bic_nl) or np.isnan(bic_bl)): 
                        bic_nl, bic_bl = compare_fits(source, survey, "H gamma", file_nl, file_bl)
                        if(np.isnan(bic_nl) or np.isnan(bic_bl)): 
                            bic_nl, bic_bl = compare_fits(source, survey, "H delta", file_nl, file_bl)
            else: 
                bic_nl, bic_bl = compare_fits(source, survey, "H beta", file_nl, file_bl)
                if(np.isnan(bic_nl) or np.isnan(bic_bl)): 
                    bic_nl, bic_bl = compare_fits(source, survey, "H gamma", file_nl, file_bl)
                    if(np.isnan(bic_nl) or np.isnan(bic_bl)): 
                        bic_nl, bic_bl = compare_fits(source, survey, "H delta", file_nl, file_bl)
                        if(np.isnan(bic_nl) or np.isnan(bic_bl)): 
                            bic_nl, bic_bl = compare_fits(source, survey, "Mg II] 2803", file_nl, file_bl)
            if(np.isnan(bic_nl) or np.isnan(bic_bl)): 
                best_fit = 'None'
            else: 
                if(bic_nl < bic_bl): 
                    delta_bic = bic_bl - bic_nl
                    if(delta_bic >= bic_threshold): 
                        best_fit = 'nl'
                else: 
                    delta_bic = bic_nl - bic_bl 
                    if(delta_bic >= bic_threshold): 
                        best_fit = 'bl'
        else: 
            best_fit = 'None'
        fit.append(best_fit)
        ## inspect for outflow 
        if(best_fit == 'nl'): 
            file, error_message = retrieve_nburstfile(source, source['survey'], 'nl_')
            if(file != 'None'): 
                hdu = fits.open(file[0])
                h3 = hdu[1].data['H3'][0][1]
                h3_err = hdu[1].data['E_H3'][0][1]
                if(h3 + h3_err < 0): 
                    if(source['continuum'] == True): 
                        outflow_mask_nl_cont[i] = True 
                    else: 
                        outflow_mask_nl_nocont[i] = True 
        elif(best_fit == 'bl'): 
            file, error_message = retrieve_nburstfile(source, source['survey'], 'bl_')
            if(file != 'None'): 
                hdu = fits.open(file[0])
                h3 = hdu[1].data['H3'][0][1]
                h3_err = hdu[1].data['E_H3'][0][1]
                if(h3 + h3_err < 0): 
                    if(source['continuum'] == True): 
                        outflow_mask_nl_cont[i] = True 
                    else: 
                        outflow_mask_nl_nocont[i] = True
    o_nl_cont = table[outflow_mask_nl_cont]
    add_outflow(o_nl_cont, survey, True, 'nl', nl_sig_limit)
    o_nl_nocont = table[outflow_mask_nl_nocont]
    add_outflow(o_nl_nocont, survey, False, 'nl', nl_sig_limit)
    o_bl_cont = table[outflow_mask_bl_cont]
    add_outflow(o_bl_cont, survey, True, 'bl', nl_sig_limit)
    o_bl_nocont = table[outflow_mask_bl_nocont]
    add_outflow(o_bl_nocont, survey, False, 'bl', nl_sig_limit)
    table['best fit'] = fit 
    ## check that outflow fit is best 
    outflow_bool = [False]*len(table)
    for i, source in enumerate(table): 
        out_file, error_message = retrieve_nburstfile(source, survey, 'out_')
        if(out_file != 'None' and table['best fit'][i] != 'None'):
            if(table['best fit'][i] == 'nl'): 
                orig_file, error_message = retrieve_nburstfile(source, survey, 'nl_')
            elif(table['best fit'][i] == 'bl'): 
                orig_file, error_message= retrieve_nburstfile(source, survey, 'bl_')
            bic_out, bic_orig = compare_fits(source, survey, "[O III] 5008", out_file, orig_file) ## z < 0.8
            if(np.isnan(bic_out) or np.isnan(bic_orig)): 
                    bic_out, bic_orig = compare_fits(source, survey, "[O II] 3729", out_file, orig_file) ## z < 1.41 
            if(~np.isnan(bic_out) and ~np.isnan(bic_orig)): 
                if(bic_out < bic_orig): 
                    delta_bic = bic_orig- bic_out
                    if(delta_bic > bic_threshold): 
                        outflow_bool[i] = True 
    table['outflow'] = outflow_bool
    return table     

def check_errors(table, survey): 
    errors = []
    for i, source in enumerate(table):
        if(table['outflow'][i] == True): 
            file, error_message = retrieve_nburstfile(source, survey, 'out_') 
        elif(table['best fit'][i] != 'None'): 
            file, error_message = retrieve_nburstfile(source, survey, f"{table['best fit'][i]}_")
        else: 
            error_message = "[ERROR] Could not determine best fit"
        errors.append(error_message)
    table['error message'] = errors 
    return table 


def run_script(table, survey): 
    continuum, no_continuum = divide_into_continuum_groups(table, survey)
    m = [True if name in continuum['name'] else False for name in table['name']]
    table['continuum'] = [False] * len(table)
    table['continuum'][m] = True 
    for source in continuum: 
        with open(f"{output_location}/{survey}_continuum.txt", 'a') as f: 
            if(survey == 'sdss'):
                f.write(f"{sdss_filepath}/{source['source file']}\n")
            else: 
                f.write(f"{desi_filepath}/{source['desi tile']}/{source['desi lastnight']}/{source['source file']}\n")
    for source in no_continuum: 
        with open(f"{output_location}/{survey}_nocontinuum.txt", 'a') as f: 
                    if(survey == 'sdss'):
                        f.write(f"{sdss_filepath}/{source['source file']}\n")
                    else: 
                        f.write(f"{desi_filepath}/{source['desi tile']}/{source['desi lastnight']}/1d/{source['source file']}\n")
    table = perform_fitting(table, survey, 150, 2)
    table = check_errors(table, survey)
    return table 
    
    

def main():
    ## clear everything from outputs and start fresh 
    os.system(f"rm {output_location}/**")

    hdu = fits.open('./nburst_fittingscript/imbh_sample_prepared.fits')
    imbh = Table(hdu[1].data)
    m = (imbh['SDSS_CLASS'] == 'QSO') | (imbh['DESI_SPECTYPE'] == 'QSO')
    input = Table(imbh[m])
    m = (imbh['survey'] == 'sdss')
    sdss_sources = imbh[m]
    m = imbh['survey'] == 'desi'
    desi_sources = imbh[m]
    sdss_table = run_script(sdss_sources, 'sdss')
    desi_table = run_script(desi_sources, 'desi')
    sdss_table['continuum'] = sdss_table['continuum'].astype(bool)
    desi_table['continuum'] = desi_table['continuum'].astype(bool)
    sdss_table['best fit'] = sdss_table['best fit'].astype(str)
    desi_table['best fit'] = desi_table['best fit'].astype(str)
    table_tot = vstack([sdss_table, desi_table])
    write_output_tbl(table_tot)

if __name__ == "__main__":
    main()
