# CMBstat project

## Next steps
  + Method for importing data (including polarisations?)
  + Method for subsetting polygonal subarea of sky.
  + Calculate area given a subset of sky (from HEALPix pixel sizes).
    + [Done in python.](http://healpy.readthedocs.io/en/latest/generated/healpy.query_polygon.html?highlight=polygon)
  + Triangulation of sphere.
    + [Test to know if a vector is inside a spherical triangle.](https://math.stackexchange.com/questions/1175362/test-to-know-if-a-vector-is-inside-a-spherical-triangle)
  + What to do with the polarisation data (Q_STOKES, U_STOKES)?
    + [Details here: Scroll to 'Polarisation Convention / Internal Convention'](http://healpix.sourceforge.net/html/intronode6.htm)
  
## Suggestions
  + Provide HEALPix projection onto 2D surface. [See here.](http://sufoo.c.ooco.jp/program/healpix.html)
  + Provide class *skywin* to hold polygonal sky window (like *owin* in [spatstat](https://cran.r-project.org/web/packages/spatstat/index.html) package):
    + Make *summary(skywin)* return area and boundary information
  + Provide class *cmbDataFrame* that also holds Nside, units, coordinate system, ordering scheme, etc. 
    + Make *summary(cmbDataFrame)* return [number of pixels](http://healpy.readthedocs.io/en/latest/healpy_pix.html#nside-npix-resolution), area, units, coordinate system and common statistics.
  + Provide options: (1) longitude and latitude in degree, (2) longitude and co-latitude in radians. [See 'lonlat' here.](http://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.pix2ang.html#healpy.pixelfunc.pix2ang)
    + | Place     | Latitude  | Colatitude  | 
      | --------- | --------- | ----------- |
      | Nth Pole  | 90&deg;   | 0&deg;      |
      | Equator   | 0&deg;    | 90&deg;     |
      | Sth Pole  | -90&deg;  | 180&deg;    |
  + Provide common intensity unit conversions (K_CMB <-> K_RJ <-> MJy/sr). [See here.](https://irsasupport.ipac.caltech.edu/index.php?/Knowledgebase/Article/View/181/20/what-are-the-intensity-units-of-the-planck-all-sky-maps-and-how-do-i-convert-between-them)
  + Provide conversion between RING and NESTED numbering schemes for *cmbDataFrame*. [See here.](http://healpy.readthedocs.io/en/latest/healpy_pix.html#conversion-between-nested-and-ring-schemes)
    > "It is in the RING scheme that Fourier transforms with spherical harmonics are easy to implement.    
    >  NESTED tree structure allows one to implement efficiently all applications involving nearest-neighbour searches, and also allows for an immediate construction of the fast Haar wavelet transform" [Source.](http://healpix.sourceforge.net/html/intronode4.htm)
      
      

## Notes on Planck maps 
  + Source: http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html
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
  + [Data analysis methods for CMB](http://iopscience.iop.org/article/10.1088/0034-4885/70/6/R02/meta)
  + [How to make maps from CMB data without losing information](http://iopscience.iop.org/article/10.1086/310631)
  + [CMB and astrophysical component maps wiki](https://wiki.cosmos.esa.int/planckpla/index.php/CMB_and_astrophysical_component_maps)
  + [Google Sky with CMB overlay](https://www.google.com.au/sky/)
  + [Free software for viewing FITS files](https://heasarc.gsfc.nasa.gov/ftools/fv/fv_download.html)
  + [HEALPix license info](http://healpix.sourceforge.net/downloads.php)
  + [Original HEALPix paper](http://cosmocoffee.info/viewtopic.php?t=64)
  + [An errata and notes on original HEALPix paper](http://blog.tiaan.com/link/2009/09/04/healpix-errata-and-additional-notes)
  + [HEALPix information page](http://healpix.sourceforge.net/)
  + [NASA HEALPix information page](http://healpix.jpl.nasa.gov/index.shtml)
  + [NASA Planck knowledge base](https://irsasupport.ipac.caltech.edu/index.php?/Knowledgebase/List/Index/20/planck)
  + [healpy (python) documentation page](http://healpy.readthedocs.io/en/latest/index.html)
  + [HEALPix C++ documentation](http://healpix.sourceforge.net/html/Healpix_cxx/index.html)
  + [HEALPix C subroutines](http://healpix.sourceforge.net/html/csub.htm)
  + [CRAN packages in R for astronomy](https://asaip.psu.edu/forums/software-forum/459833927)
    + [R bindings for Google's s2: Spherical geometry](https://cran.r-project.org/web/packages/s2/index.html) and [github repo](https://github.com/spatstat/s2) and [C++ source](https://code.google.com/archive/p/s2-geometry-library/)
    + [FITS file manipulation](http://cran.us.r-project.org/web/packages/astro/index.html)
    + [Spherical plotting](https://cran.r-project.org/web/packages/sphereplot/)
    + [FITSio](https://cran.r-project.org/web/packages/FITSio/index.html)
    + [misc astro functions (e.g. see angSep function)](https://cran.r-project.org/web/packages/astroFns/index.html)
  + [Journal of Statistical Software](https://www.jstatsoft.org)
  + [The R Journal](https://journal.r-project.org)
  
## Methods to look into more
  + > "Standard operations of numerical analysis which one might wish to execute on the sphere include convolutions with local and  global kernels, Fourier analysis with spherical harmonics and power spectrum estimation, wavelet decomposition, nearest-neighbour searches, topological analysis, including searches for extrema or zero-crossings, computing Minkowski functionals, extraction of patches and finite differencing for solving partial differential equations." [Source: HEALPix doc.](http://healpix.sourceforge.net/html/intronode3.htm)
  + [Pixel window functions.](http://healpix.jpl.nasa.gov/html/intronode14.htm)
