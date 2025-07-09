import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from astropy.io import fits


fits_file_path = 'dla_100kpc_15.fits'
#########Open the FITS file
with fits.open(fits_file_path) as hdul:
    # Access the data from the FITS file
    data = hdul[1].data  # Access the first data extension (change the index if needed)
    print(hdul[1].header)
    wave = data['WAVE    ']
    median = data['STACK_MEDIAN ']
    mean = data['STACK_MEAN ']
    

plt.plot(wave,median,'k',label='data')
plt.xlabel('Wave')
plt.ylabel('Stacked_median')
plt.xlim(1190.,1300.)
plt.legend()
plt.tight_layout()
plt.show()


