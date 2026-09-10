import argparse
import glob
import os
import re
import subprocess
import tempfile
from collections import defaultdict
import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table
from tqdm import tqdm


SDSS_FILEPATH = '/Users/f007znp/Research/processed/SDSS_BOSS'
DESI_FILEPATH = '/Users/f007znp/Research/processed/DESI'
OUTPUT_LOCATION = '/Users/f007znp/Research/IMBH/nburst_fittingscript/outputs'
FIT_OUT_SDSS = '/Users/f007znp/Research/pro/trial4/SDSS_BOSS'
FIT_OUT_DESI = '/Users/f007znp/Research/pro/trial4/DESI'
SOFTENING_PARAM_FILE = './nburst_fittingscript/sdss_softening_param.txt'
NBURSTS_PRO = '/Users/f007znp/Research/nbursts/idl/ppxf_nbursts_c.pro'
SSP_PATH = '/Users/f007znp/Research/stellar_templates/XSL/Kroupa/'

BIN_WIDTH = 25
MAX_NUM_RUNS = 3
SIG_CEILING = 850

os.makedirs(OUTPUT_LOCATION, exist_ok=True)
sp = ascii.read(SOFTENING_PARAM_FILE)


def sdss_flux_mag(f, band):
    idx = np.where(sp['filter'] == band)[0][0]
    b = sp['b'][idx]
    return -2.5/np.log10(10) * (np.arcsinh((f/1e9)/(2*b)) + np.log10(b))

def desi_flux_mag(f):
    f = np.asarray(f, dtype=float)
    mag = np.full_like(f, np.nan)
    good = f > 0
    mag[good] = 22.5 - 2.5*np.log10(f[good])
    return mag


def sdss_lrg_definition(t):
    i = sdss_flux_mag(t['SDSS_CALIBFLUX_i'], 'i')
    z = sdss_flux_mag(t['SDSS_CALIBFLUX_z'], 'z')
    r = sdss_flux_mag(t['SDSS_CALIBFLUX_r'], 'r')
    lrg_izw = (i - z > 0.7) & (i - t['WISE_w1mpro'] > 2.143*(i - z) - 0.2) & (z < 19.95) & (i > 19.9)
    lrg_riw = (r - i > 0.98) & (r - t['WISE_w1mpro'] > 2*(r - i)) & (i - z > 0.625) & (z < 19.95) & (i > 19.9)
    return lrg_izw | lrg_riw

def desi_lrg_definition(t):
    g = desi_flux_mag(t['DESI_FLUX_G'])
    r = desi_flux_mag(t['DESI_FLUX_R'])
    z = desi_flux_mag(t['DESI_FLUX_Z'])
    onea = z - t['WISE_w1mpro'] > 0.8*(r - z) - 0.6
    oneb = ((g - t['WISE_w1mpro'] > 2.6) & (g - r > 1.4)) | (r - t['WISE_w1mpro'] > 1.8)
    onec = (r - z > (z - 16.83)*0.45) & (r - z > (z - 13.80)*0.19)
    oned = r - z > 0.7
    onee = desi_flux_mag(t['desi fiberflux_z']) < 21.5
    return onea & oneb & onec & oned & onee


def compute_bic(hdu, k, target, window=30):
    idx = np.where(hdu[2].data["LINE_ID"][0] == target)[0][0]
    wavelength = hdu[2].data['WAVE'][0][idx]
    spec = hdu[1].data
    wave, flux, fit_total, error = spec['WAVE'][0], spec['FLUX'][0], spec['FIT'][0], spec['ERROR'][0]
    mask = (wave > wavelength - window) & (wave < wavelength + window)
    flux, fit_total, error = flux[mask], fit_total[mask], error[mask]
    chi2 = np.sum(((flux - fit_total) / error) ** 2)
    return chi2 + k*np.log(len(flux)), chi2


def run_gdl_script(gdl_code, workdir=None, timeout=1800, total=None, pbar_desc="Fitting"):
    full_code = f"start.pro\n .r {NBURSTS_PRO}\n{gdl_code}\nexit\n"
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pro', delete=False) as f:
        f.write(full_code)
        script_path = f.name
    process_re = re.compile(r"Processing\s+(\d+)\s*/\s*(\d+)")
    try:
        process = subprocess.Popen(["gdl", script_path], stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True, cwd=workdir, bufsize=1)
        pbar = tqdm(total=total, desc=pbar_desc, unit="obj")
        last_i = 0
        lines = []
        for line in process.stdout:
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
        full_output = "".join(lines)
        with open(f"{OUTPUT_LOCATION}/gdl_last_run.log", "a") as f:
            f.write(full_output)
        if process.returncode != 0:
            print(f"WARNING: GDL exited with code {process.returncode} for '{pbar_desc}' — check gdl_last_run.log")
        return full_output
    finally:
        os.remove(script_path)


def parse_ids(filename, survey):
    base = filename.split("/")[-1]
    if survey == 'sdss':
        # spec-PLATE-MJD-FIBERID.fits
        parts = base.split("-")
        return (parts[1].lstrip("0"), parts[2].lstrip("0"),
                parts[3].split(".")[0].lstrip("0"))
    else:
        # desiSp1d_TILE-thruNIGHT-FIBER.fits(.gz)
        prefix_and_tile, night_part, fiber_part = base.split("-")
        tile = prefix_and_tile.split("_")[1]
        night = night_part.replace("thru", "")
        fiber = fiber_part.split(".")[0]
        return (tile.lstrip("0"), night.lstrip("0"), fiber.lstrip("0"))


def build_source_path(base, survey, filename, tile=None, night=None):
    if survey == 'desi':
        return f"{base}/{tile}/{night}/1d/{filename}"
    else:
        return f"{base}/{filename}"


def find_nburst_file(base_out, stage, survey, id1, id2, id3):
    pattern = f"{base_out}/{stage}_/nbursts_{survey}_*{id1}*_*{id2}*_*{id3}*"
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f"No output file matched: {pattern}")
    return matches[0]


def bin_value(x, width=BIN_WIDTH):
    x = float(np.asarray(x).item()) if np.ndim(x) > 0 else float(x)
    return int(round(x / width) * width)


class Fitting:
    def __init__(self, filelist_path, survey, fit_stellarpop, nl_sig_limit):
        self.file = filelist_path
        self.survey = survey
        self.fit_stellarpop = fit_stellarpop
        self.nl_sig_limit = nl_sig_limit
        self.outpath = FIT_OUT_SDSS if survey == 'sdss' else FIT_OUT_DESI

    def fit_nl(self, total):
        common = (f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=2,emexcl=0,"
                  f"emlt1=[1],start=[0,100,0,80,3000,-1.2],/force_sigma,")
        if self.fit_stellarpop:
            gdl_code = (common +
                f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4,siglimits=[200,{self.nl_sig_limit}],"
                f"path_ssp='{SSP_PATH}',prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/nl_'")
        else:
            gdl_code = (common + "/disable_stpop,"
                f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4,siglimits=[200,{self.nl_sig_limit}],"
                f"outpath='{self.outpath}/nl_'")
        run_gdl_script(gdl_code, workdir="/Users/f007znp/Research/pro", total=total, pbar_desc="NL batch")

    def fit_bl(self, sig_max, total):
        common = (f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=3,emexcl=0,"
                  f"emlt1=[1],emlt2=[2],start=[0,100,0,80,0,400,3000,-1.2],/force_sigma,")
        if self.fit_stellarpop:
            gdl_code = (common +
                f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4,siglimits=[200,{self.nl_sig_limit},{sig_max}],"
                f"path_ssp='{SSP_PATH}',prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/bl_'")
        else:
            gdl_code = (common + "/disable_stpop,"
                f"lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4,siglimits=[200,{self.nl_sig_limit},{sig_max}],"
                f"outpath='{self.outpath}/bl_'")
        run_gdl_script(gdl_code, workdir="/Users/f007znp/Research/pro", total=total, pbar_desc="BL batch")

    def fit_outflow(self, classification, start_velo, sig_max_o, sig_max, total):
        addstart3 = f"addstart=[0,0,0,0,-{start_velo},0,0,0],"
        addstart4 = f"addstart=[0,0,0,0,-{start_velo},0,0,0,0,0],"

        if classification in (0, -1):
            base = (f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=3,emexcl=0,"
                    f"emlt1=[1,2],start=[0,100,0,80,0,150,3000,-1.2],/force_sigma,")
            siglim = f"siglimits=[200,150,{sig_max_o}],symlosvdid=[2],noh4h6losvdid=[2],"
            if self.fit_stellarpop:
                gdl_code = (base + addstart3 +
                    "lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4," + siglim +
                    f"path_ssp='{SSP_PATH}',prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/out_'")
            else:
                gdl_code = (base + "/disable_stpop," + addstart3 +
                    "lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4," + siglim +
                    f"outpath='{self.outpath}/out_'")
        else:
            base = (f"process_survey,'{self.survey}',inptable='{self.file}',nlosvd=4,emexcl=0,"
                    f"emlt1=[1,2],emlt2=[3],start=[0,100,0,80,0,150,0,400,3000,-1.2],/force_sigma,")
            siglim = f"siglimits=[200,150,{sig_max_o},{sig_max}],symlosvdid=[2],noh4h6losvdid=[2],"
            if self.fit_stellarpop:
                gdl_code = (base + addstart4 +
                    "lammin=3700,lammax=9000,degree=2,mdegree=5,moments=4," + siglim +
                    f"path_ssp='{SSP_PATH}',prefix='SB_',suffix='_XSL_Kroupa_PC.fits',outpath='{self.outpath}/out_'")
            else:
                gdl_code = (base + "/disable_stpop," + addstart4 +
                    "lammin=3700,lammax=9000,degree=-1,mdegree=5,moments=4," + siglim +
                    f"outpath='{self.outpath}/out_'")
        run_gdl_script(gdl_code, workdir="/Users/f007znp/Research/pro", total=total, pbar_desc="Outflow batch")



def classify_source(hdu_nl, hdu_bl, k_nl, k_bl, bl_sig_min, sig_max):
    sig = hdu_bl[1].data['SIG'][0][2]
    if sig < bl_sig_min or sig > sig_max:
        return 0, 9999, 0

    lammin, lammax = hdu_nl[0].header['LAMMIN'], hdu_nl[0].header['LAMMAX']
    fit_failed = 0
    bic_nl = bic_bl = np.nan

    if lammax > 6562.8 and lammin < 6562.8:
        target = "H alpha"
        idx = np.where(hdu_nl[2].data['LINE_ID'][0] == target)[0][0]
        nl_fl, nl_err = hdu_nl[2].data['FLUX'][0][idx], hdu_nl[2].data['FLUX_ERR'][0][idx]
        idx = np.where(hdu_bl[2].data['LINE_ID'][0] == target)[0][0]
        bl_fl, bl_err = hdu_bl[2].data['FLUX'][0][idx], hdu_bl[2].data['FLUX_ERR'][0][idx]

        if np.isnan(nl_fl) or np.isnan(bl_fl):
            fit_failed = 1
        else:
            bic_nl = compute_bic(hdu_nl, k_nl, target)[0] if nl_fl/nl_err >= 3 else np.nan
            bic_bl = compute_bic(hdu_bl, k_bl, target)[0] if bl_fl/bl_err >= 3 else np.nan
    else:
        def line_vals(hdu, target):
            idx = np.where(hdu[2].data['LINE_ID'][0] == target)[0][0]
            return hdu[2].data['FLUX'][0][idx], hdu[2].data['FLUX_ERR'][0][idx]

        nl_mg_fl, nl_mg_err = line_vals(hdu_nl, "Mg II] 2803")
        bl_mg_fl, bl_mg_err = line_vals(hdu_bl, "Mg II] 2803")
        nl_hb_fl, nl_hb_err = line_vals(hdu_nl, "H beta")
        bl_hb_fl, bl_hb_err = line_vals(hdu_bl, "H beta")

        mg_ok = not (np.isnan(nl_mg_fl) or np.isnan(bl_mg_fl))
        hb_ok = not (np.isnan(nl_hb_fl) or np.isnan(bl_hb_fl))

        bic_nl_mg = bic_bl_mg = bic_nl_hb = bic_bl_hb = np.nan
        if mg_ok:
            bic_nl_mg = compute_bic(hdu_nl, k_nl, "Mg II] 2803")[0] if nl_mg_fl/nl_mg_err >= 3 else np.nan
            bic_bl_mg = compute_bic(hdu_bl, k_bl, "Mg II] 2803")[0] if bl_mg_fl/bl_mg_err >= 3 else np.nan
        if hb_ok:
            bic_nl_hb = compute_bic(hdu_nl, k_nl, "H beta")[0] if nl_hb_fl/nl_hb_err >= 3 else np.nan
            bic_bl_hb = compute_bic(hdu_bl, k_bl, "H beta")[0] if bl_hb_fl/bl_hb_err >= 3 else np.nan

        if not mg_ok and not hb_ok:
            fit_failed = 1
        elif not mg_ok and hb_ok:
            bic_nl, bic_bl = bic_nl_hb, bic_bl_hb
        elif mg_ok and not hb_ok:
            bic_nl, bic_bl = bic_nl_mg, bic_bl_mg
        elif mg_ok and hb_ok:
            bic_nl = np.nanmean([bic_nl_mg, bic_nl_hb])
            bic_bl = np.nanmean([bic_bl_mg, bic_bl_hb])

    if fit_failed == 1 or np.isnan(bic_nl) or np.isnan(bic_bl):
        return 0, 9999, 1

    if bic_nl < bic_bl:
        delta = bic_bl - bic_nl
        classification = 0 if delta > 6 else -1
    else:
        delta = bic_nl - bic_bl
        classification = 1 if delta > 6 else -1
    return classification, delta, 0


def check_outflow(hdu_out, original_hdu, k_out, k_original):
    lammin, lammax = hdu_out[0].header['LAMMIN'], hdu_out[0].header['LAMMAX']
    if not (lammin <= 5008 and lammax >= 5008):
        return 0, 0

    target = "[O III] 5008"
    idx = np.where(original_hdu[2].data['LINE_ID'][0] == target)[0][0]
    fl, fl_err = original_hdu[2].data['FLUX'][0][idx], original_hdu[2].data['FLUX_ERR'][0][idx]
    bic_orig = compute_bic(original_hdu, k_original, target)[0] if fl/fl_err >= 3 else np.nan

    idx = np.where(hdu_out[2].data['LINE_ID'][0] == target)[0][0]
    fl, fl_err = hdu_out[2].data['FLUX'][0][idx], hdu_out[2].data['FLUX_ERR'][0][idx]
    bic_out = compute_bic(hdu_out, k_out, target)[0] if fl/fl_err >= 3 else np.nan

    if np.isnan(bic_orig) or np.isnan(bic_out):
        return 0, 1

    delta_out = abs(bic_orig - bic_out)
    return (1 if delta_out > 6 else 0), 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str)
    parser.add_argument("output", type=str)
    args = parser.parse_args()

    input_tbl = Table(fits.getdata(args.input, 1))

    if not os.path.exists(args.output) or os.path.getsize(args.output) == 0:
        with open(args.output, 'w') as f:
            f.write("file,ra,dec,redshift,survey,continuum,classification,"
                     "delta_BIC,outflow,fit_fail,hit_limit,nburst_file\n")

    has_sdss = (input_tbl['SDSS_SPECOBJID'] != '') & (input_tbl['SDSS_CLASS'] != 'STAR')
    has_desi = (input_tbl['DESI_TARGETID'] != '') & (input_tbl['DESI_SPECTYPE'] != 'STAR')
    keep = (has_sdss | has_desi) & (input_tbl['Rec_Z'] > 0)

    is_desi_highz_nonLRG = (
        (input_tbl['Rec_Z'] > 0.7) & (input_tbl['survey'] == 'desi') & (~desi_lrg_definition(input_tbl))
    )
    continuum_col = np.where(is_desi_highz_nonLRG, 0, 1)

    
    state = {}
    for i in np.where(keep)[0]:
        row = input_tbl[i]
        survey = row['survey']
        filename = row['source file']
        id1, id2, id3 = parse_ids(filename, survey)
        state[i] = {
            'file': filename, 'survey': survey,
            'ra': row['ra'], 'dec': row['dec'], 'redshift': row['Rec_Z'],
            'continuum': int(continuum_col[i]),
            'nl_sig_limit': 150.0, 'num_runs': 0,
            'finished': False, 'hit_limit': 0,
            'id1': id1, 'id2': id2, 'id3': id3,
        }

    output_lines = []
    round_num = 0

    while any(not s['finished'] for s in state.values()):
        round_num += 1
        active = {k: s for k, s in state.items() if not s['finished']}

       
        nl_groups = defaultdict(list)
        for k, s in active.items():
            nl_groups[(s['survey'], s['continuum'], s['nl_sig_limit'])].append(k)

        for (survey, cont, sig_limit), keys in nl_groups.items():
            base = SDSS_FILEPATH if survey == 'sdss' else DESI_FILEPATH
            batch_path = f"{OUTPUT_LOCATION}/batch_nl_r{round_num}_{survey}_{cont}_{int(sig_limit)}.txt"
            with open(batch_path, 'w') as f:
                for k in keys:
                    s = state[k]
                    path = build_source_path(base, survey, s['file'], s['id1'], s['id2'])
                    f.write(f"{path}\n")

            Fitting(batch_path, survey, bool(cont), sig_limit).fit_nl(total=len(keys))

            fit_out = FIT_OUT_SDSS if survey == 'sdss' else FIT_OUT_DESI
            for k in keys:
                s = state[k]
                nl_file = find_nburst_file(fit_out, 'nl', survey, s['id1'], s['id2'], s['id3'])
                hdu_nl = fits.open(nl_file)
                s['nl_file'] = nl_file
                s['nl_sig'] = hdu_nl[1].data['SIG'][0][-1]
                s['hdu_nl'] = hdu_nl
                s['k_nl'] = hdu_nl[0].header['NWLFIT'] - hdu_nl[0].header['DOF']

       
        for s in active.values():
            s['bl_sig_min'] = 1.5 * s['nl_sig']
            s['sig_max'] = 5 * s['nl_sig']
            s['sig_max_bin'] = bin_value(s['sig_max'])

        bl_groups = defaultdict(list)
        for k, s in active.items():
            bl_groups[(s['survey'], s['continuum'], s['nl_sig_limit'], s['sig_max_bin'])].append(k)

        for (survey, cont, sig_limit, sig_max_bin), keys in bl_groups.items():
            base = SDSS_FILEPATH if survey == 'sdss' else DESI_FILEPATH
            batch_path = f"{OUTPUT_LOCATION}/batch_bl_r{round_num}_{survey}_{cont}_{int(sig_limit)}_{sig_max_bin}.txt"
            with open(batch_path, 'w') as f:
                for k in keys:
                    s = state[k]
                    path = build_source_path(base, survey, s['file'], s['id1'], s['id2'])
                    f.write(f"{path}\n")

            Fitting(batch_path, survey, bool(cont), sig_limit).fit_bl(sig_max_bin, total=len(keys))

            fit_out = FIT_OUT_SDSS if survey == 'sdss' else FIT_OUT_DESI
            for k in keys:
                s = state[k]
                bl_file = find_nburst_file(fit_out, 'bl', survey, s['id1'], s['id2'], s['id3'])
                hdu_bl = fits.open(bl_file)
                s['bl_file'] = bl_file
                s['hdu_bl'] = hdu_bl
                s['k_bl'] = hdu_bl[0].header['NWLFIT'] - hdu_bl[0].header['DOF']

       
        for s in active.values():
            classification, delta, fit_failed = classify_source(
                s['hdu_nl'], s['hdu_bl'], s['k_nl'], s['k_bl'], s['bl_sig_min'], s['sig_max']
            )
            s['classification'] = classification
            s['delta'] = delta
            s['fit_failed'] = fit_failed
            s['outflow'] = 0

            if fit_failed != 1:
                if classification == 1:
                    h3 = s['hdu_bl'][1].data['H3'][0][1]
                    h3_err = s['hdu_bl'][1].data['E_H3'][0][1]
                    sig = s['hdu_bl'][1].data['SIG'][0][1]
                else:
                    h3 = s['hdu_nl'][1].data['H3'][0][1]
                    h3_err = s['hdu_nl'][1].data['E_H3'][0][1]
                    sig = s['hdu_bl'][1].data['SIG'][0][1]  # preserved from original — verify intentional
                s['h3'], s['h3_err'], s['sig_for_outflow'] = h3, h3_err, sig
                s['needs_outflow'] = (h3 + h3_err) < 0
            else:
                s['needs_outflow'] = False

        
        outflow_active = {k: s for k, s in active.items() if s['needs_outflow']}
        for s in outflow_active.values():
            s['start_velo'] = 1.5 * s['sig_for_outflow']
            s['sig_max_o'] = 3 * s['sig_for_outflow']
            s['sig_max_o_bin'] = bin_value(s['sig_max_o'])
            s['start_velo_bin'] = bin_value(s['start_velo'], width=10)

        out_groups = defaultdict(list)
        for k, s in outflow_active.items():
            out_groups[(s['survey'], s['continuum'], s['classification'],
                        s['start_velo_bin'], s['sig_max_o_bin'], s['sig_max_bin'])].append(k)

        for (survey, cont, classif, start_velo_bin, sig_max_o_bin, sig_max_bin), keys in out_groups.items():
            base = SDSS_FILEPATH if survey == 'sdss' else DESI_FILEPATH
            batch_path = (f"{OUTPUT_LOCATION}/batch_out_r{round_num}_{survey}_{cont}_{classif}_"
                          f"{start_velo_bin}_{sig_max_o_bin}_{sig_max_bin}.txt")
            with open(batch_path, 'w') as f:
                for k in keys:
                    s = state[k]
                    path = build_source_path(base, survey, s['file'], s['id1'], s['id2'])
                    f.write(f"{path}\n")

            Fitting(batch_path, survey, bool(cont), state[keys[0]]['nl_sig_limit']).fit_outflow(
                classif, start_velo_bin, sig_max_o_bin, sig_max_bin, total=len(keys)
            )

            fit_out = FIT_OUT_SDSS if survey == 'sdss' else FIT_OUT_DESI
            for k in keys:
                s = state[k]
                out_file = find_nburst_file(fit_out, 'out', survey, s['id1'], s['id2'], s['id3'])
                hdu_out = fits.open(out_file)
                s['out_file'] = out_file
                k_out = hdu_out[0].header['NWLFIT'] - hdu_out[0].header['DOF']

                if s['classification'] in (0, -1):
                    original_hdu, k_original = s['hdu_nl'], s['k_nl']
                else:
                    original_hdu, k_original = s['hdu_bl'], s['k_bl']

                outflow, of_fit_failed = check_outflow(hdu_out, original_hdu, k_out, k_original)
                s['outflow'] = outflow
                if of_fit_failed:
                    s['fit_failed'] = 1

        
        for k, s in active.items():
            if s['outflow'] == 1:
                s['nburst_file'] = s['out_file']
            elif s['classification'] == 1:
                s['nburst_file'] = s['bl_file']
            else:
                s['nburst_file'] = s['nl_file']

            hdu_final = fits.open(s['nburst_file'])
            sig_final = hdu_final[1].data['SIG'][0][1][0]

            if np.isclose(sig_final, s['nl_sig_limit'], atol=1e-1) and s['num_runs'] < MAX_NUM_RUNS:
                s['nl_sig_limit'] = min(2*s['nl_sig_limit'], SIG_CEILING)
                s['num_runs'] += 1
                s['hit_limit'] = 1
            else:
                s['finished'] = True
                output_lines.append(
                    f"{s['file']},{s['ra']},{s['dec']},{s['redshift']},{s['survey']},{s['continuum']},"
                    f"{s['classification']},{s['delta']},{s['outflow']},{s['fit_failed']},"
                    f"{s['hit_limit']},{s['nburst_file']}\n"
                )

    with open(args.output, 'a') as f:
        f.writelines(output_lines)


if __name__ == "__main__":
    main()