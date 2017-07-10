from astropy.io import fits

data, header = fits.getdata("../../CMB_map_smica1024.fits", header=True)

fits.writeto('CMB_testmap_1024_10cols.fits', data[:10], header, overwrite=True)