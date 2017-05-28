# CMBstat project

## Next steps
  + Method for importing data (including polarisations?)
  + Method for subsetting polygonal subarea of sky.
  + Calculate area given a subset of sky (from HEALPix pixel sizes).
  
## Suggestions
  + Provide HEALPix projection onto 2D surface. [See here.](http://sufoo.c.ooco.jp/program/healpix.html)
  + Provide class *skywin* to hold polygonal sky subarea:
    + Make *summary(skywin)* return area and boundary information.
  + Provide options: (1) longitude and latitude in degree, (2) longitude and co-latitude in radians. [See 'lonlat' here.](http://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.pix2ang.html#healpy.pixelfunc.pix2ang)
    + | Place     | Latitude  | Colatitude  | 
      | --------- | --------- | ----------- |
      | Nth Pole  | 90&deg;   | 0&deg;      |
      | Equator   | 0&deg;    | 90&deg;     |
      | Sth Pole  | -90&deg;  | 180&deg;    |
      
      

## Notes on plank maps 
  + Source: http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html
  + All Sky Maps are in HEALPix format, with Nside 1024 or 2048, in Galactic coordinates, and NESTED ordering. 
  + Signal given in units of Kcmb for 30-353 GHz (microwave is in this band).
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
  + [HEALPix information page](http://healpix.sourceforge.net/)
  + [healpy (python) documentation page](http://healpy.readthedocs.io/en/latest/index.html)
  + [HEALPix C++ documentation](http://healpix.sourceforge.net/html/Healpix_cxx/index.html)
  + [CRAN packages in R for astronomy](https://asaip.psu.edu/forums/software-forum/459833927)
    + [FITS file manipulation](http://cran.us.r-project.org/web/packages/astro/index.html)
    + [Spherical plotting](https://cran.r-project.org/web/packages/sphereplot/)
    + [FITSio](https://cran.r-project.org/web/packages/FITSio/index.html)
    + [misc astro functions (e.g. see angSep function)](https://cran.r-project.org/web/packages/astroFns/index.html)
  + [Journal of Statistical Software](https://www.jstatsoft.org)
  + [The R Journal](https://journal.r-project.org)
  
## Notes on analysis
  + > Standard operations of numerical analysis which one might wish to execute on the sphere include convolutions with local and  global kernels, Fourier analysis with spherical harmonics and power spectrum estimation, wavelet decomposition, nearest-neighbour searches, topological analysis, including searches for extrema or zero-crossings, computing Minkowski functionals, extraction of patches and finite differencing for solving partial differential equations. [Source: HEALPix doc.](http://healpix.sourceforge.net/html/intronode3.htm)
