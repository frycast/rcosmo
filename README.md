# CMBstat project

## Next steps
  + Method for importing data (including polarisations?)
  + Method for subsetting polygonal subarea of sky.
  + Calculate area given a subset of sky (from HEALPix pixel sizes).

## Notes on plank maps 
  + Source: http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html
  + All Sky Maps are in HEALPix format, with Nside 1024 or 2048, in Galactic coordinates, and nested ordering. 
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
  + [Journal of Statistical Software](https://www.jstatsoft.org)
  + [The R Journal](https://journal.r-project.org)
