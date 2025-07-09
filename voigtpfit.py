import Lya_zelda as Lya
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from radiative_corner import make_corner_plots
# from vpfit import lym_radiative
import VoigtFit





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


def lym_radiative(Geometry,data_1,z_avg,lo,l1,FWHM_t):

    
    your_grids_location = '/home/astro-s5/dla/DLA/Grids/'
    Lya.funcs.Data_location = your_grids_location
    
    LyaRT_Grid = Lya.load_Grid_Line( Geometry)
    
    data = pd.read_csv(f'{data_1}')
    wave = np.array(data['Wavelength']) * (1 + z_avg)
    flux = np.array(data['Flux']) / (1 + z_avg)
    err = np.array(data['Error']) / (1 + z_avg)
    # id_x = np.where((wave>1213.2201 * (1 + z_avg))&(wave< 1218.0801 * (1 + z_avg)))[0]
    id_x = np.where((wave >= lo * (1 + z_avg))&(wave<= l1 * (1 + z_avg)))[0]
    # Filter the data based on the indices
    w_Arr = wave[id_x]
    print(w_Arr / (1+z_avg))
    f_Arr = np.array(data['Flux'])[id_x] / (1 + z_avg)  # Adjust flux accordingly
    s_Arr = np.array(data['Error'])[id_x] / (1 + z_avg)  # Adjust error accordingly
    # flux_low,flux_med,flux_up = vpfit (wave,z_avg,low_nh1,med_nh1,up_nh1)
    
    ##### id_tau
    id_tau = np.where((wave > 1215.6701 * (1 + z_avg))&(wave <= l1 * (1 + z_avg) ))[0]
    # Integrate the flux over the wavelength range using the trapezoidal rule
    integrated_flux = np.trapz(flux[id_tau], wave[id_tau])  # in erg/s/cm^2

    # Defining the model parameters:
    z_t      = z_avg  # redshift of the source
    V_t      = 40.0  # Outflow expansion velocity [km/s]
    log_N_t  = 20.   # Logarithmic of the neutral hydrogen column density [cm**-2]
    t_t      = 0.01  # Dust optical depth
    log_EW_t = 1.5   # Logarithmic the intrinsic equivalent width [A]
    W_t      = 0.5   # Intrinsic width of the line [A]
    # F_t      = np.max(np.array(data['Flux']))  # Total flux of the line
    F_t      = integrated_flux # Total flux of the line
    # Defining the quality of the line profile:
    PNR_t  =15. # Signal to noise ratio of the maximum of the line.
    # FWHM_t = 0.5318057564651215  # Full width half maximum diluting the line. Mimics finite resolution. [A]
    # FWHM_t = 0.53 # Full width half maximum diluting the line. Mimics finite resolution. [A]
    # FWHM_t = 0.87132904 
    # PIX_t  = (w_Arr[1]-w_Arr[0]) # Wavelength binning of the line. [A]
    # PIX_t = 0.27
    PIX_t = 0.127
    # PIX_t = 0.27
    #print(PIX_t)
    # w_Arr , f_Arr , s_Arr = Lya.Generate_a_real_line( z_t , V_t, log_N_t, t_t, F_t, log_EW_t, W_t , PNR_t, FWHM_t, PIX_t, LyaRT_Grid, Geometry )

    # w_Arr , f_Arr , s_Arr  = Lya.Generate_a_real_line( z_t , V_t, log_N_t, t_t, F_t, log_EW_t, W_t , PNR_t, FWHM_t, PIX_t, LyaRT_Grid, Geometry )

    w_pix_Arr , f_pix_Arr = Lya.plot_a_rebinned_line( w_Arr , f_Arr , PIX_t )



    N_walkers = 300 # Number of walkers
    N_burn    = 300 # Number of steps to burn-in
    N_steps   = 500 # Number of steps to run after burning-in
    MODE = 'DNN'  # Deep Neural Network Mode
    #MODE = 'PSO'    # fast Particle Swarm Optimization

    log_V_in , log_N_in , log_t_in , log_E_in , W_in , z_in , Best = Lya.MCMC_get_region_6D( MODE , w_Arr , f_Arr , s_Arr , FWHM_t , PIX_t , LyaRT_Grid , Geometry,Geometry_Mode='Outflow' )
    # print(log_V_in,log_N_in,z_in,W_in,Best)
    sampler = Lya.MCMC_Analysis_sampler_5( w_Arr , f_Arr , s_Arr , FWHM_t ,N_walkers, N_burn, N_steps , Geometry, LyaRT_Grid , z_in=z_in , log_V_in=log_V_in , 
    log_N_in=log_N_in , log_t_in=log_t_in , log_E_in=log_E_in , W_in=W_in )

    # sampler = Lya.MCMC_Analysis_sampler_5( w_Arr  , f_Arr  , s_Arr  , FWHM_t, N_walkers, N_burn, N_steps , Geometry, LyaRT_Grid , z_in=z_in, log_V_in=[1.193561598242636,2.8987564203054068], 
    #                                       log_N_in=[18.730008832647438, 21.81905213917198], log_t_in=log_t_in , log_E_in=log_E_in , W_in=W_in )


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

    

    return z_50,V_50,t_50,W_50,log_N_50,log_E_50,F_t,FWHM_t,PIX_t,LyaRT_Grid,flat_samples


