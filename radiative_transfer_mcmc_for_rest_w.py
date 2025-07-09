import Lya_zelda as Lya
# import FLaREON as Lya
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from radiative_corner import make_corner_plots
from matplotlib import rcParams
import VoigtFit
# from astropy.cosmology import Planck18
import astropy.units as u
# Set the default font family
rcParams['font.family'] = 'DejaVu Sans'

######################### Voigt fit For med #############
def vpfit(wave,nh1_min,nh1_med,nh1_up):
    f = 0.4164  # Oscillator strength for Lyman alpha
    velocity_dispersion = 4e6  # cm/s
    gamma = 6.265e8  # s^-1 total damping constant
    wave_0 = 1215.6701#*(1+z_avg)
    
    tau_min = VoigtFit.funcs.voigt.Voigt(wl=wave, l0=wave_0, f=f, N=10**nh1_min, b=velocity_dispersion, gam=gamma, z=0)
    tau_med = VoigtFit.funcs.voigt.Voigt(wl=wave, l0=wave_0, f=f, N=10**(nh1_med), b=velocity_dispersion, gam=gamma, z=0)
    tau_up =  VoigtFit.funcs.voigt.Voigt(wl=wave, l0=wave_0, f=f, N=10**(nh1_up), b=velocity_dispersion, gam=gamma, z=0)
    
    flux_low = np.exp(-tau_min)
    flux_med = np.exp(-tau_med)
    flux_up = np.exp(-tau_up)
    return flux_low,flux_med,flux_up



your_grids_location = '/home/dharmender/New_dla/Grids/'
Lya.funcs.Data_location = your_grids_location


data_1 = 'stacked_median_with_different_data.csv'
data_2 = 'weighted_mean_stacked_different_data.csv'
data_med = 'stacked_median_data.csv'
data_wm = 'weighted_mean_stacked_data.csv'
# z_avg = 2.4567479
z_avg = 2.6367

lo = 1205.0
l1 = 1225.0

data_median = pd.read_csv(f'{data_med}')
data_wmean = pd.read_csv(f'{data_wm}')


# Geometry = input("enter the Geometry Bicone_X_Slab_Out, Galactic_Wind, Thin_Shell_Cont : ")
Geometry = 'Thin_Shell_Cont'
# LyaRT_Grid = Lya.load_Grid_Line( Geometry)
LyaRT_Grid = Lya.load_Grid_Line( 'Thin_Shell_Cont' )

###----- Data_1 #######    
### Load data for analysis
data = pd.read_csv(data_1)
wave = np.array(data['Wavelength'])
flux = np.array(data['Flux'])
err = np.array(data['Error'])

low_nh1, med_nh1, up_nh1 = 21.0, 21.21, 21.5

# Filtering based on wavelength range
flux_low, flux_med, flux_up = vpfit(data_median['Wavelength'], low_nh1, med_nh1, up_nh1)

# Filtering based on wavelength range
id_x = np.where((wave >= lo) & (wave < l1))[0]
w_Arr = wave[id_x]*(1+z_avg)
f_Arr = flux[id_x]/(1+z_avg)
s_Arr = err[id_x]/(1+z_avg)
# Integrate flux in selected range
id_tau = np.where((wave >= 1215.6701) & (wave < l1))[0]
integrated_flux = np.trapz(flux[id_tau], wave[id_tau])

##-- radiative ---##
PNR = 100000.
######---------------------------------------------------------#############  

# Defining the model parameters:
z_t      = z_avg  # redshift of the source
V_t      = 60.0  # Outflow expansion velocity [km/s]
log_N_t  = 20.5   # Logarithmic of the neutral hydrogen column density [cm**-2]
t_t      = 0.01  # Dust optical depth
log_EW_t = 0.1   # Logarithmic the intrinsic equivalent width [A]
W_t      = 0.7   # Intrinsic width of the line [A]
F_t      = integrated_flux #* (1+z_avg)# Total flux of the line
# Defining the quality of the line profile:
PNR_t  =15. # Signal to noise ratio of the maximum of the line.
# FWHM_t = 0.531805756 # Full width half maximum diluting the line. Mimics finite resolution. [A]
# FWHM_t = 0.53 # Full width half maximum diluting the line. Mimics finite resolution. [A]
FWHM_t = 0.87132904 
# PIX_t  = (w_Arr[1]-w_Arr[0]) # Wavelength binning of the line. [A]
# PIX_t = 0.27 # Wavelength binning of the line. [A]
PIX_t = 0.01

w_pix_Arr , f_pix_Arr = Lya.plot_a_rebinned_line( w_Arr, f_Arr , PIX_t )
N_walkers = 200 # Number of walkers
N_burn    = 200 # Number of steps to burn-in
N_steps   = 300 # Number of steps to run after burning-in
# N_walkers = 300 # Number of walkers
# N_burn    = 300 # Number of steps to burn-in
# N_steps   = 500 # Number of steps to run after burning-in
MODE = 'DNN'  # Deep Neural Network Mode
#MODE = 'PSO'    # fast Particle Swarm Optimization

log_V_in , log_N_in , log_t_in , log_E_in , W_in , z_in , Best = Lya.MCMC_get_region_6D( MODE , w_Arr, f_Arr , s_Arr , FWHM_t , PIX_t , LyaRT_Grid , Geometry,Geometry_Mode='Outflow' )
# print(log_V_in,log_N_in,z_in,log_t_in,log_E_in,W_in,Best)
# sampler = Lya.MCMC_Analysis_sampler_5( w_Arr, f_Arr , s_Arr , FWHM_t,N_walkers, N_burn, N_steps , Geometry, LyaRT_Grid , z_in=z_in , log_V_in=log_V_in , 
#                                         log_N_in=log_N_in , log_t_in=log_t_in , log_E_in=log_E_in , W_in=W_in )


sampler = Lya.MCMC_Analysis_sampler_5( w_Arr  , f_Arr  , s_Arr  , FWHM_t, N_walkers, N_burn, N_steps , Geometry, LyaRT_Grid , z_in=z_in, log_V_in=log_V_in, 
                                      log_N_in=[21.001083265, 21.4053198], log_t_in=log_t_in , log_E_in=[-0.290366441644, 1.03702940] , W_in=W_in )

# sampler = Lya.MCMC_Analysis_sampler_5( w_Arr  , f_Arr  , s_Arr  , FWHM_t, N_walkers, N_burn, N_steps , Geometry, LyaRT_Grid , z_in=z_in, log_V_in=[1.0,3.2], 
#                                       log_N_in=[20.401083265, 21.3053198],log_t_in=[-0.01,2.01] , log_E_in=[-0.290366441644, 1.03702940] , W_in=[0.01,5.0] )


Q_Arr = [ 16 , 50 , 84 ] # You can add more percentiles here, like 95

perc_matrix_sol , flat_samples = Lya.get_solutions_from_sampler( sampler , N_walkers , N_burn , N_steps , Q_Arr )

# redshift.
z_16     =     perc_matrix_sol[ 3 , 0 ] # corresponds to Q_Arr[0]
z_50     =     perc_matrix_sol[ 3 , 1 ] # corresponds to Q_Arr[1]
z_84     =     perc_matrix_sol[ 3 , 2 ] # corresponds to Q_Arr[2]

# Expansion velocity.
V_16     = 10**perc_matrix_sol[ 0 , 0 ]
V_50     = 10**perc_matrix_sol[ 0 , 1 ]
V_84     = 10**perc_matrix_sol[ 0 , 2 ]

# dust optical depth.
t_16     = 10**perc_matrix_sol[ 2 , 0 ]
t_50     = 10**perc_matrix_sol[ 2 , 1 ]
t_84     = 10**perc_matrix_sol[ 2 , 2 ]

# Intrinsic width.
W_16     =     perc_matrix_sol[ 5 , 0 ]
W_50     =     perc_matrix_sol[ 5 , 1 ]
W_84     =     perc_matrix_sol[ 5 , 2 ]

# Logarithmic of the intrinsic equivalent width.
log_E_16 =     perc_matrix_sol[ 4 , 0 ]
log_E_50 =     perc_matrix_sol[ 4 , 1 ]
log_E_84 =     perc_matrix_sol[ 4 , 2 ]

# Logarithmic of the HI column density.
log_N_16 =     perc_matrix_sol[ 1 , 0 ]
log_N_50 =     perc_matrix_sol[ 1 , 1 ]
log_N_84 =     perc_matrix_sol[ 1 , 2 ]

print( 'The true redshift                 is' , z_t      , 'and the predicted is' , z_50     , '(-' , z_50-z_16         , ', +' , z_84-z_50         , ')' )
print( 'The true expansion velocity       is' , V_t      , 'and the predicted is' , V_50     , '(-' , V_50-V_16         , ', +' , V_84-V_50         , ')' )
print( 'The true dust optical depth       is' , t_t      , 'and the predicted is' , t_50     , '(-' , t_50-t_16         , ', +' , t_84-t_50         , ')' )
print( 'The true intrinsic width          is' , W_t      , 'and the predicted is' , W_50     , '(-' , W_50-W_16         , ', +' , W_84-W_50         , ')' )
print( 'The true log of HI column density is' , log_N_t  , 'and the predicted is' , log_N_50 , '(-' , log_N_50-log_N_16 , ', +' , log_N_84-log_N_50 , ')' )
print( 'The true log of equivalent width  is' , log_EW_t , 'and the predicted is' , log_E_50 , '(-' , log_E_50-log_E_16 , ', +' , log_E_84-log_E_50 , ')' )




# #### Compute line For data_1-----###########   

w_One_Arr , f_One_Arr , _  = Lya.Generate_a_real_line( z_50, V_50, log_N_50, t_50, F_t, log_E_50, W_50, PNR, FWHM_t, PIX_t, LyaRT_Grid, Geometry )
#### Make cooler
w_pix_One_Arr , f_pix_One_Arr = Lya.plot_a_rebinned_line( w_One_Arr , f_One_Arr, PIX_t )


###------- Data_2----######### 

#### Compute line
data_mean = pd.read_csv(f'{data_2}')
wave_mean = np.array(data_mean['Wavelength']) 
flux_mean = np.array(data_mean['Flux']) 
err_mean = np.array(data_mean['Error'])

# flux_mean = lum_mean / (4*np.pi* (lum_dist_cm**2))
# print("the error is :",err,err_mean)

id_x_1 = np.where((wave_mean >= lo)&(wave_mean< l1 ))[0]
# id_x = np.where((tau_mean >= 10.))[0]

# Filter the data based on the indices
w_Arr_mean = wave_mean[id_x_1] * (1 + z_avg)
# print(w_Arr / (1+z_avg))
f_Arr_mean = flux_mean[id_x_1] / (1 + z_avg)  # Adjust flux accordingly
s_Arr_mean = err_mean[id_x_1] / (1 + z_avg)  # Adjust error accordingly
# flux_low,flux_med,flux_up = vpfit (wave,z_avg,low_nh1,med_nh1,up_nh1)

##### id_tau
# id_tau_mean = np.where((wave_mean > 1215.6701)&(tau_mean>=10. ))[0]
id_tau_mean = np.where((wave_mean > 1215.6701)&(wave_mean < l1))[0]
# Integrate the flux over the wavelength range using the trapezoidal rule
integrated_flux_mean = np.trapz(flux_mean[id_tau_mean], wave_mean[id_tau_mean])  # in erg/s/cm^2

# Defining the model parameters:
z_t_mean      = z_avg  # redshift of the source
V_t_mean      = 90.0  # Outflow expansion velocity [km/s]
log_N_t_mean  = 20.5   # Logarithmic of the neutral hydrogen column density [cm**-2]
t_t_mean      = 0.01  # Dust optical depth
log_EW_t_mean = 0.1   # Logarithmic the intrinsic equivalent width [A]
W_t_mean      = 0.7   # Intrinsic width of the line [A]
F_t_mean      = integrated_flux_mean #*(1+z_avg) # Total flux of the line
print(F_t_mean)
# Defining the quality of the line profile:
PNR_t  =15. # Signal to noise ratio of the maximum of the line.
FWHM_t_mean = 0.8718057  # Full width half maximum diluting the line. Mimics finite resolution. [A]
PIX_t_mean = 0.01

w_pix_Arr_mean , f_pix_Arr_mean = Lya.plot_a_rebinned_line( w_Arr_mean, f_Arr_mean , PIX_t_mean )
##log_v_in 
log_V_in_mean , log_N_in_mean , log_t_in_mean , log_E_in_mean , W_in_mean , z_in_mean , Best_mean = Lya.MCMC_get_region_6D( MODE , w_Arr_mean , f_Arr_mean , s_Arr_mean , FWHM_t_mean , PIX_t_mean , LyaRT_Grid , Geometry)#,Geometry_Mode='Outflow' )
# print(log_V_in,log_N_in,z_in,W_in,Brest)
# sampler_mean = Lya.MCMC_Analysis_sampler_5(w_Arr_mean, f_Arr_mean , s_Arr_mean , FWHM_t_mean ,N_walkers, N_burn, N_steps , Geometry, LyaRT_Grid , z_in=z_in_mean , log_V_in=log_V_in_mean , 
#                                             log_N_in=[21.03000883, 21.339052] , log_t_in=log_t_in_mean , log_E_in=log_E_in_mean , W_in=W_in_mean )

sampler_mean = Lya.MCMC_Analysis_sampler_5(w_Arr_mean , f_Arr_mean , s_Arr_mean , FWHM_t_mean ,N_walkers, N_burn, N_steps , Geometry, LyaRT_Grid , z_in=z_in_mean , log_V_in=[1.0, 3.027929] , 
                                            log_N_in=[20.430832, 21.35198], log_t_in=log_t_in_mean , log_E_in=[-0.290366, 1.03702] , W_in=[0.1,5.05])

# sampler_mean = Lya.MCMC_Analysis_sampler_5( w_Arr_mean  , f_Arr_mean  , s_Arr_mean, FWHM_t_mean, N_walkers, N_burn, N_steps , Geometry, LyaRT_Grid , z_in=[2.634665,2.6386665], log_V_in=[0.42448994721, 3.72792402909], 
#                                        log_N_in=[21.030832, 21.505198], log_t_in=[-4.7124233348, 2.313688041008]  , log_E_in=[-0.190366441644, 0.73702940] , W_in=[ 0.02656853, 6.38552421] )


Q_Arr = [ 16 , 50 , 84 ] # You can add more percentiles here, like 95

perc_matrix_sol_mean , flat_samples_mean = Lya.get_solutions_from_sampler( sampler_mean , N_walkers , N_burn , N_steps,Q_Arr )

# # redshift.
z_16_mean     =     perc_matrix_sol_mean[ 3 , 0 ] # corresponds to Q_Arr[0]
z_50_mean     =     perc_matrix_sol_mean[ 3 , 1 ] # corresponds to Q_Arr[1]
z_84_mean     =     perc_matrix_sol_mean[ 3 , 2 ] # corresponds to Q_Arr[2]

# Expansion velocity.
V_16_mean     = 10**perc_matrix_sol_mean[ 0 , 0 ]
V_50_mean     = 10**perc_matrix_sol_mean[ 0 , 1 ]
V_84_mean     = 10**perc_matrix_sol_mean[ 0 , 2 ]

# dust optical depth.
t_16_mean     = 10**perc_matrix_sol_mean[ 2 , 0 ]
t_50_mean     = 10**perc_matrix_sol_mean[ 2 , 1 ]
t_84_mean     = 10**perc_matrix_sol_mean[ 2 , 2 ]

# Intrinsic width.
W_16_mean     =     perc_matrix_sol_mean[ 5 , 0 ]
W_50_mean     =     perc_matrix_sol_mean[ 5 , 1 ]
W_84_mean     =     perc_matrix_sol_mean[ 5 , 2 ]

# Logarithmic of the intrinsic equivalent width.
log_E_16_mean =     perc_matrix_sol_mean[ 4 , 0 ]
log_E_50_mean =     perc_matrix_sol_mean[ 4 , 1 ]
log_E_84_mean =     perc_matrix_sol_mean[ 4 , 2 ]

# Logarithmic of the HI column density.
log_N_16_mean =     perc_matrix_sol_mean[ 1 , 0 ]
log_N_50_mean =     perc_matrix_sol_mean[ 1 , 1 ]
log_N_84_mean =     perc_matrix_sol_mean[ 1 , 2 ]
# redshift.
# z_16_mean     =     perc_matrix_sol_mean[ 3  ] # corresponds to Q_Arr[0]
# z_50_mean     =     perc_matrix_sol_mean[ 3 ] # corresponds to Q_Arr[1]
# z_84_mean     =     perc_matrix_sol_mean[ 3 ] # corresponds to Q_Arr[2]

# # Expansion velocity.
# V_16_mean     = 10**perc_matrix_sol_mean[ 0 ]
# V_50_mean     = 10**perc_matrix_sol_mean[ 0  ]
# V_84_mean     = 10**perc_matrix_sol_mean[ 0  ]

# # dust optical depth.
# t_16_mean     = 10**perc_matrix_sol_mean[ 2 ]
# t_50_mean     = 10**perc_matrix_sol_mean[ 2 ]
# t_84_mean     = 10**perc_matrix_sol_mean[ 2  ]

# # Intrinsic width.
# W_16_mean     =     perc_matrix_sol_mean[ 5 ]
# W_50_mean     =     perc_matrix_sol_mean[ 5 ]
# W_84_mean     =     perc_matrix_sol_mean[ 5 ]

# # Logarithmic of the intrinsic equivalent width.
# log_E_16_mean =     perc_matrix_sol_mean[ 4 ]
# log_E_50_mean =     perc_matrix_sol_mean[ 4 ]
# log_E_84_mean =     perc_matrix_sol_mean[ 4]

# # Logarithmic of the HI column density.
# log_N_16_mean =     perc_matrix_sol_mean[ 1 ]
# log_N_50_mean =     perc_matrix_sol_mean[ 1 ]
# log_N_84_mean =     perc_matrix_sol_mean[ 1 ]

print( 'The true redshift                 is' , z_t_mean      , 'and the predicted is' , z_50_mean     , '(-' , z_50_mean-z_16_mean         , ', +' , z_84_mean-z_50_mean         , ')' )
print( 'The true expansion velocity       is' , V_t_mean      , 'and the predicted is' , V_50_mean     , '(-' , V_50_mean-V_16_mean         , ', +' , V_84_mean-V_50_mean         , ')' )
print( 'The true dust optical depth       is' , t_t_mean      , 'and the predicted is' , t_50_mean     , '(-' , t_50_mean-t_16_mean         , ', +' , t_84_mean-t_50_mean         , ')' )
print( 'The true intrinsic width          is' , W_t_mean      , 'and the predicted is' , W_50_mean     , '(-' , W_50_mean-W_16_mean         , ', +' , W_84_mean-W_50_mean         , ')' )
print( 'The true log of HI column density is' , log_N_t_mean  , 'and the predicted is' , log_N_50_mean , '(-' , log_N_50_mean-log_N_16_mean , ', +' , log_N_84_mean-log_N_50_mean , ')' )
print( 'The true log of equivalent width  is' , log_EW_t_mean , 'and the predicted is' , log_E_50_mean , '(-' , log_E_50_mean-log_E_16_mean , ', +' , log_E_84_mean-log_E_50_mean , ')' )




w_One_Arr_mean , f_One_Arr_mean , _  = Lya.Generate_a_real_line( z_50_mean, V_50_mean, log_N_50_mean, t_50_mean, F_t_mean, log_E_50_mean, W_50_mean, PNR, FWHM_t_mean, PIX_t_mean, LyaRT_Grid, Geometry )

#### Make cooler
w_pix_One_Arr_mean , f_pix_One_Arr_mean = Lya.plot_a_rebinned_line( w_One_Arr_mean , f_One_Arr_mean, PIX_t_mean )


#######-------- Data_1 ------------####
data_name = 'median'
make_corner_plots( flat_samples )
plt.tight_layout()
plt.savefig(f"{data_name}.pdf")
# plt.savefig("fitting_results_mean.pdf")
plt.show()
#######-------- Data_2 ------------####
data_name_mean = 'mean'
make_corner_plots( flat_samples_mean )
plt.tight_layout()
plt.savefig(f"{data_name_mean}.pdf")
# plt.savefig("fitting_results_mean.pdf")
plt.show()

### Plotting start from here
import matplotlib as mp
# Set font and tick parameters
mp.rcParams['font.family'] = 'serif'
mp.rcParams['xtick.major.size'] = 7
mp.rcParams['xtick.major.width'] = 1
mp.rcParams['xtick.minor.size'] = 4
mp.rcParams['xtick.minor.width'] = 1
mp.rcParams['ytick.major.size'] = 7
mp.rcParams['ytick.major.width'] = 1
mp.rcParams['ytick.minor.size'] = 4
mp.rcParams['ytick.minor.width'] = 1
mp.rcParams['axes.linewidth'] = 1
mp.rcParams['xtick.labelsize'] = 24
mp.rcParams['ytick.labelsize'] = 24
mp.rcParams['xtick.direction'] = 'in'
mp.rcParams['ytick.direction'] = 'in'



# Plotting the original data
fig, axs = plt.subplots(1, 2, figsize=(30, 12))

# Plotting the original data from 1210 to 1220 Å For data_1
id_w = np.where((w_pix_One_Arr>= lo*(1+z_avg)) &(w_pix_One_Arr<l1*(1+z_avg)))[0]
# id_w = np.where((tau>=10.))[0]
print(w_pix_One_Arr[id_w])
axs[0].step(data_median['Wavelength'], data_median['Flux'], where='mid', color='k', lw=2.50, label='DLA (log(N$_{HI}$/cm$^{-2})>21.0$)')
axs[0].errorbar(data_median['Wavelength'], data_median['Flux'], yerr=data_median['Error'],fmt='none', ecolor='grey', capsize=6)

####--- Overlay the MCMC results
axs[0].plot(w_pix_One_Arr[id_w]/(1+z_avg), f_pix_One_Arr[id_w]*(1+z_avg), color='r', lw=2.60, label='Ly$\\alpha$ Model')

# axs[0].plot(w_pix_One_Arr[id_w] / (1 + z_avg), f_pix_One_Arr[id_w] * (1 + z_avg), color='r', lw=2.60, label='Ly$\\alpha$ Model')

###-------- VoigtFit plot--------########
axs[0].plot(data_median['Wavelength'], (flux_low * 6.9), color='g', linestyle='--', lw=1.8)
# axs[0].plot(data_median['Wavelength'], (flux_med * 6.9), color='b', linestyle='--', lw=1.8)
# axs[0].plot(data_median['Wavelength'], (flux_up * 6.9), color='g', linestyle='--', lw=1.8)
axs[0].hlines(y=0.0, xmin=(np.min(wave)), xmax=(np.max(wave)), 
              color='k', linestyle='--', lw=1.8)
# axs[0].hlines(y=0.0, xmin=(np.min(wave / (1 + z_avg))), xmax=(np.max(wave / (1 + z_avg))), 
#               color='k', linestyle='--', lw=1.8)

# Add hatched region
hatch_start,hatch_end = 1213.4835,1217.8019
hatch_bottom,hatch_top = (np.min(flux)) - 0.070,(np.min(flux)) - 0.035

axs[0].fill_between([hatch_start, hatch_end], hatch_bottom, hatch_top, 
                    color='none', edgecolor='g', hatch='//', linewidth=0)

axs[0].set_xlabel("DLA's rest-frame wavelength (${\\AA}$)", fontsize=28)
axs[0].set_ylabel('F$_\\lambda$ (10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)',fontsize=28) #fontsize=mp.rcParams['ytick.labelsize'])
# axs[0].set_ylabel('L$_\\lambda$ (10$^{40}$ erg s$^{-1} \\AA^{-1}$)')

# Adding text to the hatched region
text_x = (hatch_start + hatch_end) / 2
text_y = (hatch_bottom + hatch_top) / 2
axs[0].text(text_x, text_y, '$\\tau > 10$', color='k', fontsize=mp.rcParams['ytick.labelsize'], ha='center', va='center')
# axs[0].text(0.25,-0.55, 'Median', transform=plt.gca().transAxes, ha='center', fontsize=mp.rcParams['ytick.labelsize'])
axs[0].set_title('Median', y=0.97, pad=-14,fontsize=27)
axs[0].set_xlim(1210., 1220.)
axs[0].set_ylim((np.min(flux)) - 0.06, 0.25)

# axs[0].set_ylim((np.min(flux) * (1 + z_avg)) - 0.06, 0.25)
# axs[0].set_title('Median')
axs[0].legend(loc=3,fontsize=28)#fontsize=mp.rcParams['ytick.labelsize']

#### -------- Data_2 ---------##############
id_w_1 = np.where((w_pix_One_Arr_mean>= lo*(1+z_avg) ) &(w_pix_One_Arr_mean<l1*(1+z_avg)))[0]
print(w_pix_One_Arr_mean[id_w_1])
# id_w_1 = np.where((w_pix_One_Arr_mean>=lo*(1+z_avg))&(w_pix_One_Arr_mean < l1*(1+z_avg)))[0]


axs[1].step(data_wmean['Wavelength'], data_wmean['Flux'], where='mid', color='k', lw=2.50, label='DLA (log(N$_{HI}$/cm$^{-2})>21.0$)')

axs[1].errorbar(data_wmean['Wavelength'], data_wmean['Flux'], yerr=data_wmean['Error'], fmt='none', ecolor='grey', capsize=6)

# Overlay the MCMC results
axs[1].plot(w_pix_One_Arr_mean[id_w_1]/(1+z_avg), f_pix_One_Arr_mean[id_w_1]*(1+z_avg), color='r', lw=2.60, label='Ly$\\alpha$ Model')

# axs[1].plot(w_pix_One_Arr_mean[id_w_1] / (1 + z_avg), f_pix_One_Arr_mean[id_w_1] * (1 + z_avg), color='r', lw=2.60, label='Ly$\\alpha$ Model')
axs[1].plot(data_wmean['Wavelength'], (flux_low * 6.9), color='g', linestyle='--', lw=1.8)
# axs[1].plot(data_median['Wavelength'], (flux_med * 6.9), color='b', linestyle='--', lw=1.8)
# axs[1].plot(data_median['Wavelength'], (flux_up * 6.9), color='g', linestyle='--', lw=1.8)
axs[1].hlines(y=0.0, xmin=(np.min(wave_mean)), xmax=(np.max(wave_mean)), color='k', linestyle='--', lw=1.8)

# axs[1].hlines(y=0.0, xmin=(np.min(wave_mean / (1 + z_avg))), xmax=(np.max(wave_mean / (1 + z_avg))), color='k', linestyle='--', lw=1.8)

# Add hatched region
hatch_bottom_mean,hatch_top_mean = (np.min(flux_mean)) - 0.070,(np.min(flux_mean)) - 0.035

axs[1].fill_between([hatch_start, hatch_end], hatch_bottom_mean, hatch_top_mean, 
                    color='none', edgecolor='g', hatch='//', linewidth=0)

# Adding text to the hatched region
text_y_mean = (hatch_bottom_mean + hatch_top_mean) / 2
axs[1].text(text_x, text_y_mean, '$\\tau > 10$', color='k', fontsize=mp.rcParams['ytick.labelsize'], ha='center', va='center')
# Add text at the top
axs[1].tick_params(axis='y', left=False, labelleft=False, right=False, labelright=False)
axs[1].spines['left'].set_visible(True)
axs[1].set_title('3$\\sigma$-clipped weighted mean', y=0.97, pad=-16,fontsize=25)
axs[1].set_xlabel("DLA's rest-frame wavelength (${\\AA}$)", fontsize=28)
axs[1].set_xlim(1210., 1220.)
axs[1].set_ylim((np.min(flux_mean)) - 0.06, 0.25)
axs[1].legend(loc=3,fontsize=mp.rcParams['ytick.labelsize'])

# Reduce the space between the plots
plt.subplots_adjust(wspace=0.03)

plt.tight_layout()
plt.savefig('spectrum_fit_mcmc_median_mean.pdf')
plt.show()









