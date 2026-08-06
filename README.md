# nburst_fittingscript

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
  device,retain=2,decomposed=0 <code>
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