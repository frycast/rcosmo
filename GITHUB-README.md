# The `rcosmo` project

### The vignette (version 1) is available by [clicking here](rcosmoVignette.pdf)
### Documentation (version 0.1, needs work) is available by [clicking here](rcosmo.pdf). 


The Cosmic Microwave Background (CMB) is remnant electromagnetic radiation from the epoch of recombination. It is the most ancient important source of data about the early universe and the key to unlocking the mysteries of the Big Bang and the structure of time and space. Spurred on by a wealth of satellite data, intensive investigations in the past few years have resulted in many physical and mathematical results to characterize CMB radiation. An advanced R programming toolkit is needed to help statisticians perform CMB data analytics. 

`rcosmo` addresses various data processing and statistical analysis needs for the present generation of CMB experiments. These needs fall into the following broad categories:
+ Importing and transforming HEALPix data in convenient CMBDataFrames
+ Geometric tools
+ Statistical tools
+ Visualisation

The current version of `rcosmo` includes the following functionality:
+	Generation of a comprehensive data frame of CMB observations, which include HEALPix indices, metadata, CMB intensities and their
  corresponding spherical and/or cartesian coordinates, as well as polarization data
+	Window subsetting tools for investigating circular, convex and non-convex polygonal sub-regions on the sphere
+	Fast empirical covariance and variogram estimation
+	Implementation and analysis of spherical harmonics, spherical wavelets, etc
+	Various methods for CMB map visualization, such as interactive 3D full sky maps rendered with OpenGL, Mollweide projection and HEALPix boundary plotting
+ Spherical geometry tools such as shortest distance between two points, calculate spherical angles, shortest distance between a point and a region, etc


## Installation (Windows and OSX)
First install the devtools package:
```
install.packages("devtools")
```
Then, if using Microsoft Windows, install the latest RTools from CRAN at this link [here](https://cran.r-project.org/bin/windows/Rtools/) if it is not already installed on your device.
Then use devtools to install rcosmo:
```
devtools::install_github("VidaliLama/rcosmo")
```
If you did not install RTools and are using RStudio then you will be prompted to install RTools. After installing RTools you should run `devtools::install_github("VidaliLama/rcosmo")` again.

## Installation (Linux)

```
install.packages("devtools")
library(devtools)
find_rtools()
source("https://raw.githubusercontent.com/r-lib/remotes/master/install-github.R")$value("r-lib/remotes")
remotes::install_github("VidaliLama/rcosmo")
```

## Next steps
  + Kriging (e.g. for equatorial region)
  + 2D projection [Mollweide view](https://en.wikipedia.org/wiki/Mollweide_projection)
  
## Suggestions
  + Provide common intensity unit conversions (K_CMB <-> K_RJ <-> MJy/sr). [See here.](https://irsasupport.ipac.caltech.edu/index.php?/Knowledgebase/Article/View/181/20/what-are-the-intensity-units-of-the-planck-all-sky-maps-and-how-do-i-convert-between-them)

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

#### Data
  + [Planck legacy archive](http://pla.esac.esa.int/pla/#home)
#### Related papers and wikis
  + [FITS Standard](https://fits.gsfc.nasa.gov/fits_standard.html)
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
#### Journals
  + [Journal of Statistical Software](https://www.jstatsoft.org)
  + [The R Journal](https://journal.r-project.org)
#### Source files for stars and constellations
  + [https://www.iau.org/public/themes/constellations/](https://www.iau.org/public/themes/constellations/)
  + [http://www.midnightkite.com/index.aspx?AID=0&URL=StarChartFAQ](http://www.midnightkite.com/index.aspx?AID=0&URL=StarChartFAQ)
  + [http://pbarbier.com/constellations/boundaries.html](http://pbarbier.com/constellations/boundaries.html)
#### Other examples of HEALPix data
  + [http://healpix.sourceforge.net/resources.php](http://healpix.sourceforge.net/resources.php)
  + [http://aladin.u-strasbg.fr/java/Demo/AladinDemo.gml](http://aladin.u-strasbg.fr/java/Demo/AladinDemo.gml)
  + [http://www.atnf.csiro.au/people/mcalabre/WCS/example_data.html](http://www.atnf.csiro.au/people/mcalabre/WCS/example_data.html)

