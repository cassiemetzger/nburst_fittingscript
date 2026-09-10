# nburst_fittingscript
Running this script requires knowledge of Nbursts, GDL, IDLAstro, and several stellar templates. Please see "Installing dependencies" to learn how to use Nbursts. 

This spectral fitting code proceeds as follows: 
<p align="center">
<img width="600" height="500" alt="Spectral_fitting_flowchart" src="plots/spectral_fitting_simple.png" />
</p>
Spectra are grouped by redshift. All spectrum under redshift 0.7 are fit with a stellar continuum </li>
For spectra above 0.7, we investigate whether or not they are an LRG. SDSS sources are classified as LRGs according to the criteria defined in <a href = "https://arxiv.org/abs/1508.04473"> Dawson+2016</a>. DESI sources are classified as LRGs according to the criteria outlined in <a href = "https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..181Z/abstract"> Zhou+2020</a>.

A narrow line + broad line fit is then applied. A 60 Angstrom window around Halpha is defined and the BIC is computed for the NL + BL fit and then again for NL + BL fit with the BL component subtracted. If Halpha falls outside the observed wavelength range, [Mg II]2796 and Hbeta are inspected instead. If $\Delta BIC > 6$, we take the fit with the lower BIC to be the best fit. Otherwise, we assume only a narrow line profile.  

Once the spectrum has been classified as either broad or narrow, the H3 component of the Gauss-Hermite function is inspected. If it is negative, we refit the spectrum with an added outflow component. We then compare the outflow fit to the previous best fit by inspecting the BIC of both fits. Again, if $\Delta BIC > 6$, we take the difference between the fits to be significant and declare the fit with the lower BIC to be the true best fit. If $\Delta BIC \leq 6$, we assume that the outflow component is extraneous and discard it. 

## Installing dependencies 
To run this code, Nbursts must be installed via Bitbucket. An <a href="https://www.atlassian.com/try/cloud/signup?bundle=bitbucket">Atlassian</a> account is required to do this. If this is your first time using Bitbucket, remember to set up an API token on your local machine! 

While in the Nbursts directory, run <code> git switch autofit </code> to ensure you're on the proper and most updated branch. 

Once Nbursts is installed, <a href = "https://gal-04.voxastro.org/~chil/Data/NBursts_models/template/">stellar population templates </a> must be downloaded. For this code, you'll need: 
<ul> 
    <li><code>pegase.hr/newgrid2_*</code></li>
    <li><code> pegase.hr/MILES/</code></li>
    <li><code> Vazdekis</code></li>
    <li><code> XSL</code></li>
</ul>
You'll want to place these templates in some folder along your Nbursts path. For instance, both Nbursts and a folder called <code>stellar_templates</code> are at the same level on my machine. 

Now, you'll want to install <code>IDLAstro</code>. 
Run <code>clone git://github.com/wlandsman/IDLAstro.git</code>. 
I placed this folder at the same level as Nbursts and stellar_templates. 

Finally, we need to set up GDL.

If you don't have it already, install <a href="https://www.macports.org/install.php">MacPorts</a>. 
Run <code>sudo port selfupdate</code> to ensure everything is up-to-date (and test that the installation was successful). 
Then, run <code>sudo port install gnudatalanguage </code>. 
Type <code>gdl</code> to test if GDL is operational. 

Then, to connect GDL and Nbursts, create a directory in your home called <code>.idl</code>. Within that directory, create a file called <code>start.pro</code>. Within <code>start.pro</code>, write the following: 

<code>on_error,2
  !path=!path+':'+EXPAND_PATH('+/Users/f007znp/Research/idlastro/')+$
              ':'+EXPAND_PATH('+/Users/f007znp/Research/nbursts/')
  device,retain=2,decomposed=0 </code>
  
Change the paths for <code>IDLAstro</code> and Nbursts as necessary. 

At the same level as Nbursts and stellar_templates, create a directory called "pro" and one called "processed". This is where your fitting inputs and outputs will go. 

## Run a test fit 
To confirm that everything is working, let's run a test. 

Navigate into your /pro/ directory and run <code>gdl start.pro</code>. You should see a number of files being compiled. Once this is completed, try opening a graphics window: <code>window,0,xs=1500,ys=600 & loadct,39</code>. 

Now, we'll download one SDSS spectrum as practice. We'll choose a file with <code>PLATE=1751 MJD=53377 FIBERID=147</code> and perform a three component fit (stellar continuum, narrow line, broad line). 

<code>get_process_sdss,1751l,53377l,147l,/plot,nlosvd=3,emexcl=0,emlt1=[1],emlt2=[2],$ start=[0,100,0,80,0,400,3000,-1.2],/force_sigma,lammin=3700,lammax=9000,$ degree=2,mdegree=5,$ path_ssp='/Users/f007znp/Research/stellar_templates/XSL/Kroupa/',$ prefix='SB_',suffix='_XSL_Kroupa_PC.fits'</code>

Let's break this call down: 

<code>/plot</code>: Plots the fit 

<code>nlosvd</code>: The number of components in the fit (in this case, 3). 

<code>emexcl</code>: The component assigned to the stellar continuum

<code>emlt1</code>: The component assigned to the narrow line

<code>emlt2</code>: The component assigned to the broad line

<code>start</code>: The start vector. It is ordered as follows: [starting velocity of component 1, starting sigma of component 1, velo of component 2, sigma of component 2, velo of component 3, sigma of component 3, age, metallicity]. The velocity and velocity dispersion are described in km/s. The age is in Myrs. 

<code>force_sigma</code>: This forces Nbursts to use the sigmas you set in your start vector. Without it, default values are used. 

<code>lammin, lammax</code>: Descibe the wavelength range over which you want to apply the fit 

<code>degree</code>: The degree of the additive Legendre polynomial used to correct the template continuum. You can set degree=-1 to remove the stellar continuum

<code>mdegree</code>: The degree of the multiplicative Legendre polynomial.

The rest of the call tells nbursts what stellar template to use and where to find it. 

This spectrum will automatically be downloaded into processed/SDSS_BOSS. The output file with the fit information will be stored in pro/SDSS_BOSS/results. 

## Fit multiple files as once
To fit multiple SDSS files at one time, you'll need to create a file that contains the location of the spectrum and the plate number of the spectrum, separated by a tab. 
For example, the first few lines might look like: 

<code> /Users/f007znp/Research/processed/SDSS_BOSS/spec-015059-59199-4555598157.fits 15059
/Users/f007znp/Research/processed/SDSS_BOSS/spec-7095-56625-0722.fits 7095
/Users/f007znp/Research/processed/SDSS_BOSS/spec-015239-59293-6036976581.fits 15239
</code>

Next, you'll want to add the survey's name and <code>inptable</code> to your call as follows: 

<code>process_survey,'sdss',inptable='../IMBH/sdss_files.txt',/plot,nlosvd=3,emexcl=0,emlt1=[1],emlt2=[2],start=[0,100,0,80,0,400,3000,-1.2],/force_sigma,lammin=3700,lammax=9000,degree=2,mdegree=5,path_ssp='/Users/f007znp/Research/stellar_templates/XSL/Kroupa/',prefix='SB_',suffix='_XSL_Kroupa_PC.fits'</code>

### Other useful details 
A narrow line only fit will look like this: <code>nlosvd=2,emexcl=0,emlt1=[1],start=[0,100,0,80,3000,-1.2]</code>

A narrow line + outflow fit will look like this: <code>nlosvd=3,emexcl=0,emlt1=[1,2],start=[0,100,0,80,0,150,3000,-1.2]</code>

A narrow line, broad line, and outflow fit will look like this: <code>nlosvd=4,emexcl=0,emlt1=[1,2], emlt2=[2], start=[0,100,0,80,0,150,0,400,3000,-1.2]</code>

To fit Gauss-Hermite polynomials, you can add in the keyword <code>moments</code>. Set <code>moments=4</code> to fit H3, H4 and set <code>moments=6</code> to fit H3,H4,H5,H6. 

You can specify the output path of the file with <code>outpath=''</code>

You can fit only a specific range of lines in your input text file with <code> imin = , imax= ,</code>. For example, <code>imin=1193,imax=1193</code> will fit only line 1193 of the input file. 

## Running the code 
You'll want to begin by running <code>prepare_file.py</code> to retrieve all necessary SDSS and DESI information. To run this file, you'll need to have <a href = "https://www.sdss4.org/dr17/spectro/spectro_access/"><code>specObj-dr17.fits</code></a>, <a href="https://sdss.org/dr19/data_access/get_data/"><code>spAll-v6_1_3.fits</code>, and <a href = "https://data.desi.lbl.gov/doc/organization/"><code>zall-tilecumulative-iron.fits</code></a>. <b>Please change the paths in the file to match the location of these files on your machine.</b> Then, run <code>prepare_file.py</code> by entering
<code> python prepare_file.py {YOUR INPUT FILE}.fits {YOUR OUTPUT FILE}.fits </code>

Next, to retrieve SDSS data, run <code>sdss_download.py</code>. You can do this by entering <code>python sdss_download.py {YOUR PREPARED FILE}.fits</code>. <b>Be sure to edit the <code>SDSS_DESTINATION</code> parameter so your data can be located</b>. 

Given the download time of DESI data, you'll need to retrieve that on your own (sorry) :/ 

Now, you're ready to run <code>nburst.scripy.py</code>. To do so, enter <code>python nburst_script.py {YOUR PREPARED INPUT FILE} output.txt</code>.  