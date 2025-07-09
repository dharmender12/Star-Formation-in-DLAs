import Lya_zelda as Lya
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18
import astropy.units as u

from luminosity_function import sfr
from radiative_transfer_mcmc import radiative_transfer

z_avg = 2.656
# l0 = float(input("Enter Tau>10 l0:"))
# l1 = float(input("Enter Tau>10 l1:"))


# l0 = 1213.4835
# l1 = 1217.8019
# z_avg = 2.526
# data = pd.read_csv('integrated_data.csv')
# data_1= pd.read_csv('stacked_median_data.csv')
# data_2 = pd.read_csv('weighted_mean_stacked_data.csv')
data_1= pd.read_csv('stacked_median_with_different_data.csv')
data_2 = pd.read_csv('weighted_mean_stacked_different_data.csv')
low_nh1 = np.median(data_1['NHI_low'])
med_nh1 = np.median(data_1['NHI_best'])
up_nh1 = np.median(data_1['NHI_up'])
# Calculate median SFR and luminosity
sfr_med, lum_med = sfr(data_1)
# Save to text file with names
with open('median_sfr_data.txt', 'w') as f:
    f.write(f'sfr_med: {sfr_med}\n')
    f.write(f'lum_med: {lum_med}\n')

# Calculate mean SFR and luminosity
sfr_mean, lum_mean = sfr(data_2)
# Save to text file with names
with open('mean_sfr_data.txt', 'w') as f:
    f.write(f'sfr_mean: {sfr_mean}\n')
    f.write(f'luminosity_mean: {lum_mean}\n')
# sfr_med,lum_med = sfr(data_1)
# np.savetxt('median_sfr_data.txt',[sfr_med,lum_med])
# sfr_mean,lum_mean = sfr(data_2)
# np.savetxt('mean_sfr_data.txt',[sfr_mean,lum_mean])

# radiative_transfer(data_1,data_2,z_avg,l0,l1,low_nh1,med_nh1,up_nh1)

# Geometry = input("enter the Geometry Bicone_X_Slab_Out, Galactic_Wind, Thin_Shell_Cont : ")
# data_median = 'stacked_median_with_different_data.csv'
# data_mean = 'weighted_mean_stacked_different_data.csv'


