# rcosmo project


### Warning: The functions in Rcosmo are still under development.
### Expected version 1 completion date is 01/07/2018. 

The Cosmic Microwave Background (CMB) is remnant electromagnetic radiation from the epoch of recombination. It is the most ancient important source of data about the early universe and the key to unlocking the mysteries of the Big Bang and the structure of time and space. Spurred on by a wealth of satellite data, intensive investigations in the past few years have resulted in many physical and mathematical results to characterize CMB radiation. An advanced R programming toolkit is needed to help statisticians perform CMB data analytics. 

`rcosmo` addresses various data processing and statistical analysis needs for the present generation of CMB experiments. These needs fall into the following broad categories:
+ Importing and transforming HEALPix data in convenient CMBDataFrames
+ Geometric tools
+ Statistical tools
+ Visualisation

The current version of rcosmo covers includes the following functionality:
+	Generation of a comprehensive data frame of CMB observations, which include HEALPix indices, metadata, CMB intensities and their
  corresponding spherical and/or cartesian coordinates, as well as polarization data
+	Window subsetting tools for investigating circular, convex and non-convex polygonal sub-regions on the sphere
+	Fast empirical covariance and variogram estimation
+	Implementation and analysis of spherical harmonics, spherical wavelets, etc
+	Various methods for CMB map visualization, such as interactive 3D full sky maps rendered with OpenGL, Mollweide projection and HEALPix boundary plotting
+ Other basic statistical and geometrical tools

## Installation
First install the devtools package:
```
install.packages("devtools")
```
Then install the latest RTools from CRAN at this link [here](https://cran.r-project.org/bin/windows/Rtools/) if it is not already installed.
Then use devtools to install rcosmo:
```
devtools::github_install("VidaliLama/rcosmo")
```
If you did not install RTools and are using RStudio then you will be prompted to install RTools. After installing RTools you should run `devtools::github_install("VidaliLama/rcosmo")` again.

## Next steps
  + Investigate R packages [geosphere](https://cran.r-project.org/web/packages/geosphere/index.html) and [VecStatGraph3D](https://www.rdocumentation.org/packages/VecStatGraphs3D/versions/1.6)
  + What to do with the polarisation data (Q_STOKES, U_STOKES).
    + [Details here: Scroll to 'Polarisation Convention / Internal Convention'](http://healpix.sourceforge.net/html/intronode6.htm)
  + Sample mean, variance and higher moments (need normalising factor with area)
  + Sample covariance and variogram (assume Gaussianity and homogeneity, don't assume isotropy):
  
    <a href="https://www.codecogs.com/eqnedit.php?latex=\frac{1}{\&hash;\text{terms&space;in&space;sum}}\sum_{t_1&space;\in&space;S^2}\sum_{\{d(t_1,t_2)&space;=&space;r\}}&space;X(t_1)X(t_2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{1}{\&hash;\text{terms&space;in&space;sum}}\sum_{t_1&space;\in&space;S^2}\sum_{\{d(t_1,t_2)&space;=&space;r\}}&space;X(t_1)X(t_2)" title="\frac{1}{\#\text{terms in sum}}\sum_{t_1 \in S^2}\sum_{\{d(t_1,t_2) = r\}} X(t_1)X(t_2)" /></a>
    
  + Kriging (e.g. for equatorial region with milky way interference)
  + 2D projection [Mollweide view](https://en.wikipedia.org/wiki/Mollweide_projection)
  
## Suggestions
  + Provide 2D HEALPix projection. [See here](http://sufoo.c.ooco.jp/program/healpix.html) and Figure 5 from [this paper](cosmocoffee.info/arxivref.php?abs=astro-ph/0409513).
  + Look at [spatstat](https://cran.r-project.org/web/packages/spatstat/index.html) package):
  + Make *summary(CMBDataFrame)* return [number of pixels](http://healpy.readthedocs.io/en/latest/healpy_pix.html#nside-npix-resolution), area, units, coordinate system and common statistics.
  + Provide common intensity unit conversions (K_CMB <-> K_RJ <-> MJy/sr). [See here.](https://irsasupport.ipac.caltech.edu/index.php?/Knowledgebase/Article/View/181/20/what-are-the-intensity-units-of-the-planck-all-sky-maps-and-how-do-i-convert-between-them)
  + Provide conversion between RING and NESTED numbering schemes for *CMBDataFrame*. [See here.](http://healpy.readthedocs.io/en/latest/healpy_pix.html#conversion-between-nested-and-ring-schemes)
    > "It is in the RING scheme that Fourier transforms with spherical harmonics are easy to implement.    
    >  NESTED tree structure allows one to implement efficiently all applications involving nearest-neighbour searches, and also allows for an immediate construction of the fast Haar wavelet transform" [Source.](http://healpix.sourceforge.net/html/intronode4.htm)
      
      

## Notes on Planck maps 
  + [official plank colourmap](https://github.com/zonca/paperplots/raw/master/data/Planck_Parchment_RGB.txt)
  + [style guide plots](https://github.com/zonca/paperplots)
  + [Source of maps](http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html).
  + All Sky Maps are in [HEALPix](http://healpix.sourceforge.net/html/intro.htm) format, with [Nside](http://healpix.sourceforge.net/html/intronode4.htm) 1024 or 2048, in Galactic coordinates, and [NESTED](http://healpix.sourceforge.net/html/intronode4.htm) ordering. [Source.](http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/)
  + Signal given in units of [Kcmb](https://irsasupport.ipac.caltech.edu/index.php?/Knowledgebase/Article/View/181/20/what-are-the-intensity-units-of-the-planck-all-sky-maps-and-how-do-i-convert-between-them) for 30-353 GHz (microwave is in this band).
  + Unpolarized maps have 2 planes: I_Stokes (intensity) and TMASK.
  + Polarized maps have 5 planes: I_Stokes (intensity), Q_Stokes and U_Stokes (linear polarization), PMASK and TMASK.
  + File names look like:
    + COM_CMB_IQU-smica_1024_R2.02_full.fits
    + COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits
    >  'R2.02' indicates Nside 1024, at 10 arcmin resolution, with polarisation.      
    >  'R2.01' indicates Nside 2048, at 5 arcmin resolution, with intensity only.       
    >  'SMICA' indicates SMICA pipeline (others: COMMANDER, NILC, SEVEM).
    
## Useful links
#### Related papers and wikis
  + [Data analysis methods for CMB](http://iopscience.iop.org/article/10.1088/0034-4885/70/6/R02/meta)
  + [How to make maps from CMB data without losing information](http://iopscience.iop.org/article/10.1086/310631)
  + [CMB and astrophysical component maps wiki](https://wiki.cosmos.esa.int/planckpla/index.php/CMB_and_astrophysical_component_maps)
  + [NASA Planck knowledge base](https://irsasupport.ipac.caltech.edu/index.php?/Knowledgebase/List/Index/20/planck)
#### Related software and apps
  + [Google Sky with CMB overlay](https://www.google.com.au/sky/)
  + [Free software for viewing FITS files](https://heasarc.gsfc.nasa.gov/ftools/fv/fv_download.html)
#### HEALPix
  + [Original paper](https://arxiv.org/pdf/astro-ph/0409513.pdf) and [discussion](http://cosmocoffee.info/viewtopic.php?t=64).
  + [An important errata and notes on original paper](http://blog.tiaan.com/link/2009/09/04/healpix-errata-and-additional-notes)
  + [License info](http://healpix.sourceforge.net/downloads.php)
  + [Information page](http://healpix.sourceforge.net/)
  + [NASA information page](http://healpix.jpl.nasa.gov/index.shtml)
  + [healpy (python) documentation page](http://healpy.readthedocs.io/en/latest/index.html)
  + [C++ documentation](http://healpix.sourceforge.net/html/Healpix_cxx/index.html)
  + [Installing HEALPix (NASA)](https://healpix.jpl.nasa.gov/html/install.htm)
  + [C subroutines](http://healpix.sourceforge.net/html/csub.htm)
#### R packages for spherical/atro data
  + See more details on all the following packages [here](https://github.com/VidaliLama/cmbstat/blob/master/PackagesForSphericalDataAnalytics.pdf)
    + [R bindings for Google's s2: Spherical geometry](https://cran.r-project.org/web/packages/s2/index.html) and [github repo](https://github.com/spatstat/s2) and [C++ source](https://code.google.com/archive/p/s2-geometry-library/)
    + [SpherWave](http://dasan.sejong.ac.kr/~dhkim/main/research/pub/SpherWaveR.pdf)
    + [CircNNTSR](https://cran.r-project.org/web/packages/CircNNTSR/index.html)
    + [sphereplot](https://cran.r-project.org/web/packages/sphereplot/)
    + [VecStatGraphs3D](https://www.rdocumentation.org/packages/VecStatGraphs3D/versions/1.6)
    + [CRAN packages in R for astronomy](https://asaip.psu.edu/forums/software-forum/459833927)
      + [astro](http://cran.us.r-project.org/web/packages/astro/index.html).
      + [FITSio](https://cran.r-project.org/web/packages/FITSio/index.html).
      + [astroFns](https://cran.r-project.org/web/packages/astroFns/index.html).
      + more...
    + [sm](https://cran.r-project.org/web/packages/sm/index.html)
    + [Directional](https://cran.r-project.org/web/packages/Directional/index.html)
    + [SphericalCubature](https://cran.r-project.org/web/packages/SphericalCubature/index.html)
    + [geosphere](https://cran.r-project.org/web/packages/geosphere/index.html)
    + [circular](https://cran.r-project.org/web/packages/circular/index.html)
#### Creating R packages and using Rcpp
  + [S4 Generics in 15 Pages, More or Less](https://www.stat.auckland.ac.nz/S-Workshop/Gentleman/S4Objects.pdf)
  + [Rcpp, Advanced R, book by Hadley Wickham](http://adv-r.had.co.nz/Rcpp.html#rcpp-package)
  + [Making R packages, book by Hadley Wickham](http://r-pkgs.had.co.nz/intro.html)
  + [Rcpp Gallery](http://gallery.rcpp.org/)
  + [Adding Rcpp to a package that uses roxygen2](https://ironholds.org/blog/adding-rcpp-to-an-existing-r-package-documented-with-roxygen2/)
  + [Installing package from private GitHub repo using PAT](https://github.com/jennybc/happy-git-with-r/blob/master/81_github-api-tokens.Rmd)
  + [R startup, Renviron details, Efficient R book](https://csgillespie.github.io/efficientR/r-startup.html)
  + [Great R packages tutorial which includes python and C++](http://www.mjdenny.com/R_Package_Pictorial.html)
#### Python related packages
  + [Astropy, editing FITS](http://www.astropy.org/astropy-tutorials/FITS-header.html) and [more](http://www.astropy.org/).
  + [rpython, run python code from R](http://rpython.r-forge.r-project.org/)
  + [rpy2-R in Python](https://rpy2.bitbucket.io/)
  + [Jupyter with the IR Kernel: Using R with Jupyter Notebooks](http://blog.revolutionanalytics.com/2015/09/using-r-with-jupyter-notebooks.html)
#### Journals for submission
  + [Journal of Statistical Software](https://www.jstatsoft.org)
  + [The R Journal](https://journal.r-project.org)
  
## Vague suggestions
  + > "Standard operations of numerical analysis which one might wish to execute on the sphere include convolutions with local and  global kernels, Fourier analysis with spherical harmonics and power spectrum estimation, wavelet decomposition, nearest-neighbour searches, topological analysis, including searches for extrema or zero-crossings, computing Minkowski functionals, extraction of patches and finite differencing for solving partial differential equations." [Source: HEALPix doc.](http://healpix.sourceforge.net/html/intronode3.htm)
  + [Pixel window functions.](http://healpix.jpl.nasa.gov/html/intronode14.htm)
