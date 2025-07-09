import numpy as np
import Lya_zelda as Lya
import pandas as pd
import matplotlib.pyplot as plt

your_grids_location = '/home/dharmender/New_dla/Grids/'
Lya.funcs.Data_location = your_grids_location

# Geometry = 'Thin_Shell_Cont'
# Geometry = 'Galactic_Wind'

# Geometry = input("enter the Geometry Bicone_X_Slab_Out, Galactic_Wind, Thin_Shell or Thin_Shell_Cont: ")


Geometry= 'Thin_Shell_Cont'

print('The Geometry used is :', Geometry)

LyaRT_Grid = Lya.load_Grid_Line( Geometry )



z_avg =0.00001#2.63678
lo = 1205.0
l1= 1225.0


data_1 = 'stacked_median_with_different_data.csv'
data_2 = 'weighted_mean_stacked_different_data.csv'
data_med = 'stacked_median_data.csv'
data_wm = 'weighted_mean_stacked_data.csv'


data = pd.read_csv(f'{data_1}')

wave = np.array(data['Wavelength']) 
flux = np.array(data['Flux']) 
# flux = lum / (4*np.pi* (lum_dist_cm**2))
err = np.array(data['Error']) 
# PIX_t = w_Arr[1] - w_Arr[0]
PIX_t = 0.01
FWHM_t = 0.87

# F_t = np.sum(f_Arr)
F_t = 0.12
# print(F_t)
# F_t = 1.0 # flux in the line
id_x = np.where((wave > lo )&(wave< l1))[0]
w_Arr = wave[id_x] 
# print(w_Arr)
f_Arr = flux[id_x]   # Adjust flux accordingly
s_Arr = err[id_x]  # Adjust error accordingly



w_pix_Arr , f_pix_Arr = Lya.plot_a_rebinned_line( w_Arr , f_Arr , PIX_t )


plt.plot( w_pix_Arr  , f_pix_Arr  )
plt.xlabel('wavelength[A]' , size=15 )
plt.ylabel('Flux' , size=15 )

# plt.ylabel('Flux density [a.u.]' , size=15 )
plt.xlim(1210,1220.)
plt.show()
machine_data =  Lya.Load_NN_model( 'Outflow' )

machine    = machine_data['Machine' ]
w_rest_Arr = machine_data[ 'w_rest' ]

# print(machine,w_rest_Arr)
### ----------- Using the DNN in the un-perturbed line profile -------------########## 
print('Using the DNN in the un-perturbed line profile')

Sol , z_sol = Lya.NN_measure( w_Arr , f_Arr , s_Arr , FWHM_t , PIX_t , machine , w_rest_Arr , N_iter=None )

print( 'The measured redshift                                                     is' , z_sol    )
print( 'The measured logarithm of the expasion velocity                           is' , Sol[0,1] )
print( 'The measured logarithm of the HI column density                           is' , Sol[0,2] )
print( 'The measured logarithm of the dust optical depth                          is' , Sol[0,3] )
print( 'The measured logarithm of the intrinsic equivalent width                  is' , Sol[0,4] )
print( 'The measured logarithm of the intrinsic            width                  is' , Sol[0,5] )
print( 'The measured shift of the true Lya wavelgnth from the maximum of the line is' , Sol[0,0] )

# PNR = 100000. # let's put infinite signal to noise in the model line

# PNR = np.max(f_Arr) / np.max(s_Arr)
PNR = 10000.0
V_sol    = 10**Sol[0,1] # Expansion velocity km/s
logN_sol =     Sol[0,2] # log of HI column density cm**-2
t_sol    = 10**Sol[0,3] # dust optical depth
logE_sol =     Sol[0,4] # log intrinsic EW [A]
W_sol    = 10**Sol[0,5] # intrinsic width [A]

# creates the line

w_One_Arr , f_One_Arr , _  = Lya.Generate_a_real_line( z_sol , V_sol, logN_sol, t_sol, F_t, logE_sol, W_sol, PNR, FWHM_t, PIX_t, LyaRT_Grid, Geometry )

# plot the target and the predicted line

w_pix_One_Arr , f_pix_One_Arr = Lya.plot_a_rebinned_line( w_One_Arr , f_One_Arr , PIX_t )

plt.plot( w_pix_Arr     , f_pix_Arr     , label='Target' )
plt.plot( w_pix_One_Arr , f_pix_One_Arr , label='1 iter' )

plt.legend(loc=0)
plt.xlabel('wavelength[A]' , size=15 )
plt.ylabel('Flux' , size=15 )
plt.xlim(1205,1222)
plt.show()


######-------------- Using the DNN with Monte Carlo perturbations -----#########


print('Using the DNN with Monte Carlo perturbations ')
Sol , z_sol , log_V_Arr , log_N_Arr , log_t_Arr , z_Arr , log_E_Arr , log_W_Arr = Lya.NN_measure( w_Arr , f_Arr , s_Arr , FWHM_t , PIX_t , machine , w_rest_Arr , N_iter=5000 )

# Redshitft
z_50     = np.percentile(    z_Arr , 50 )
z_16     = np.percentile(    z_Arr , 16 )
z_84     = np.percentile(    z_Arr , 84 )

# Expansion velocity
V_50     = 10 ** np.percentile( log_V_Arr , 50 )
V_16     = 10 ** np.percentile( log_V_Arr , 16 )
V_84     = 10 ** np.percentile( log_V_Arr , 84 )

# Logarithmic of HI column density
log_N_50 = np.percentile( log_N_Arr , 50 )
log_N_16 = np.percentile( log_N_Arr , 16 )
log_N_84 = np.percentile( log_N_Arr , 84 )

# Dust optical depth
t_50     = 10 ** np.percentile( log_t_Arr , 50 )
t_16     = 10 ** np.percentile( log_t_Arr , 16 )
t_84     = 10 ** np.percentile( log_t_Arr , 84 )

# Logarithmic of intrinsic equivalent width
log_E_50 = np.percentile( log_E_Arr , 50 )
log_E_16 = np.percentile( log_E_Arr , 16 )
log_E_84 = np.percentile( log_E_Arr , 84 )

# Intrinsic width
W_50     = 10 ** np.percentile( log_W_Arr , 50 )
W_16     = 10 ** np.percentile( log_W_Arr , 16 )
W_84     = 10 ** np.percentile( log_W_Arr , 84 )



# Compute the 100 iterations line profile
w_50th_Arr , f_50th_Arr , _  = Lya.Generate_a_real_line( z_50 , V_50, log_N_50, t_50, F_t, log_E_50, W_50, PNR, FWHM_t, PIX_t, LyaRT_Grid, Geometry )

# Get cooler profiles
w_pix_50th_Arr , f_pix_50th_Arr = Lya.plot_a_rebinned_line( w_50th_Arr , f_50th_Arr , PIX_t )

# Plot
plt.plot( w_pix_Arr    , f_pix_Arr     , label='Target'   )
plt.plot( w_pix_One_Arr , f_pix_One_Arr , label='1 iter'   )
plt.plot( w_pix_50th_Arr, f_pix_50th_Arr ,c='g', label='1000 iter')

plt.legend(loc=0)
plt.xlabel('wavelength[A]' , size=15 )
plt.ylabel('Flux' , size=15 )
plt.xlim(1205,1222)
plt.show()

print( 'The predicted redshift                 is',  z_50     , '(-' , z_50-z_16         , ', +' , z_84-z_50         , ')' )
print( 'The predicted expansion velocity       is' , V_50     , '(-' , V_50-V_16         , ', +' , V_84-V_50         , ')' )
print( 'The predicted dust optical depth       is' , t_50     , '(-' , t_50-t_16         , ', +' , t_84-t_50         , ')' )
print( 'The predicted intrinsic width          is' , W_50     , '(-' , W_50-W_16         , ', +' , W_84-W_50         , ')' )
print( 'The predicted log of HI column density is' , log_N_50 , '(-' , log_N_50-log_N_16 , ', +' , log_N_84-log_N_50 , ')' )
print( 'The predicted log of equivalent width  is' , log_E_50 , '(-' , log_E_50-log_E_16 , ', +' , log_E_84-log_E_50 , ')' )
