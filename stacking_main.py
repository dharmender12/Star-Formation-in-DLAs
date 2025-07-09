import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.stats import bootstrap
from specstack import regrid, renorm, stack
import VoigtFit
from scipy.optimize import curve_fit
from matplotlib import rcParams
from astropy.cosmology import Planck18
import astropy.units as u

# Set the default font family
rcParams['font.family'] = 'DejaVu Sans'
# Load the DLA data
# dla_data_with_cnr = pd.read_csv('filtered_dla_data_original.csv')
# dla_data_with_cnr = pd.read_csv('filtered_dla_data_cnr_4_nh1_20.8.csv')

# dla_data_with_cnr = pd.read_csv('coadd_spec_nh1_21.csv')
# dla_data_with_cnr = pd.read_csv('coadd_spec_nh1_208_21.csv')
dla_data_with_cnr = pd.read_csv('checked_merged_file.csv')

# z_dla = dla_data_with_cnr['zCNN']
# z_metal = dla_data_with_cnr['z_metal']
# z_choice = dla_data_with_cnr['z_choice']
lo = 1215.6701 - 20. 
l1 = 1215.6701 + 20.
# wave_rest_fixed = np.arange(lo,l1,0.26)
# wave_rest_fixed = np.arange(lo,l1,0.2699)
wave_rest_fixed = np.arange(lo,l1,0.2699)
# wave_rest_fixed = np.arange(1215.6701 - 0.50, 1215.6701 +0.50, 0.2699)
sigma = 3.
value_1, value_2 = 21.0,5.0
# value_1, value_2 = 20.8,4.0

# filtered_data = dla_data_with_cnr[(dla_data_with_cnr['logNHI'] >= value_1) & (dla_data_with_cnr['CNR'] >= value_2)]
filtered_data = dla_data_with_cnr[(dla_data_with_cnr['NHIfit_1_large'] >= value_1)]

print(filtered_data.columns)
z_dla = filtered_data['zCNN_1_large']

# z_metal = filtered_data['z_metal']

######################### Voigt fit #############
f = 0.4164  # Oscillator strength for Lyman alpha
velocity_dispersion = 4e6  # cm/s
gamma = 6.265e8  # s^-1 total damping constant
wave_0 = 1215.6701

# Calculate the 16th percentile lower and upper bounds for logNHI
lower_bound = np.percentile(filtered_data['NHIfit_1_large'], 16)
# lower_bound = 21.
upper_bound = np.percentile(filtered_data['NHIfit_1_large'], 84)
best_nh1 = np.percentile(filtered_data['NHIfit_1_large'], 50)#np.median(filtered_data['logNHI'])
# min_nh1 = np.min(filtered_data['logNHI'])
# min_nh1 = min(filtered_data['logNHI'].sort_values
print('min,med and max of logNHI is :',lower_bound,best_nh1,upper_bound)

# Calculate tau for Voigt profile
tau_min = VoigtFit.funcs.voigt.Voigt(wl=wave_rest_fixed, l0=wave_0, f=f, N=10**lower_bound, b=velocity_dispersion, gam=gamma, z=0)
tau_best= VoigtFit.funcs.voigt.Voigt(wl=wave_rest_fixed, l0=wave_0, f=f, N=10**best_nh1, b=velocity_dispersion, gam=gamma, z=0)
tau_max = VoigtFit.funcs.voigt.Voigt(wl=wave_rest_fixed, l0=wave_0, f=f, N=10**upper_bound, b=velocity_dispersion, gam=gamma, z=0)
flux_rest_min = np.exp(-tau_min)
flux_rest_max = np.exp(-tau_max)
flux_rest_med = np.exp(-tau_best)



print(len(filtered_data))
print("the average z_dla is : ",np.mean(z_dla))#,np.mean(z_metal))
# filepath = '/home/astro-s5/dla/DLA/Coadd_Spec/'
filepath = '/home/dharmender/New_dla/Coadd_Spec/'
# Initialize counter
counter = 0
# '/home/astro-s5/dla/DLA/Coadd_Spec'
# Collect spectra data
spectra_list = []
for i_file in range(len(filtered_data)):
    # filename = f"{filepath}{filtered_data['Filename'].iloc[i_file]}"
    # Define filename
    if filtered_data['coadd_bossname_large'].iloc[i_file] is not None and not pd.isnull(filtered_data['coadd_bossname_large'].iloc[i_file]):
        filename = f"{filepath}{filtered_data['coadd_bossname_large'].iloc[i_file]}"
    else:
        filename = f"{filepath}{filtered_data['coadd_sdssname_large'].iloc[i_file]}"
     # Use z_abs if available, otherwise use zCNN_large
    z_used = filtered_data['z_abs'].iloc[i_file] if filtered_data['z_abs'].iloc[i_file] != 0 else filtered_data['zCNN_1_large'].iloc[i_file]
    lum_dist = Planck18.luminosity_distance(z_used)
    # Convert luminosity distance to cm (1 Mpc = 3.086e24 cm)
    lum_dist_cm = lum_dist.to(u.cm).value
    # z_used = z_dla.iloc[i_file] 
    with fits.open(filename) as hdul:
        data = hdul[1].data
        wave = data['WAVE'][0]
        flux = data['FLUX'][0]
        err = data['ERROR'][0]
        conti_tm = data['CONTI_TM'][0]
        # norm_flux = flux / np.median(flux)
    
    
    try:
        # luminosity = 4 * np.pi * (lum_dist_cm ** 2) * flux *(1e-17) # in erg/s
        # luminosity =  luminosity * 10**(-40)
        wave_rest, flux_norm = renorm(wave, flux,conti_tm, z_used, l0=1215.6701 - 20., l1=1215.6701 + 20., verbose=True)
        wave_rest, flux_norm = renorm(wave,flux,conti_tm, z_used, l0=lo, l1=l1, verbose=True)
        # spectra_list.append([wave_rest, flux_norm])
        # wave_rest, luminosity_norm = renorm(wave,luminosity,conti_tm, z_used, l0=lo, l1=l1, verbose=True)
        spectra_list.append([wave_rest, flux_norm])
    except ValueError as e:
        print(f"Skipping file {filename} due to error: {e}")

# Rebin and stack spectra
rebinned, clipped = regrid(spectra_list, wave_rest_fixed,err, tau_best, sigma=sigma)
stacked_mean, std_mean, errmean, N_mean, stacked_median, std_median, ermedian, N_median,err_norm = stack(rebinned=rebinned, sigma=sigma,flux_err=err)

# print("The Error is :",err_norm)
# Define Gaussian function for curve fitting
def gaussian(x, mean, stddev):
    amplitude = 1 / (stddev * np.sqrt(2 * np.pi))
    return amplitude * np.exp(-0.5 * ((x - mean) / stddev) ** 2)

# Set font and tick parameters
plt.rcParams.update({
    'font.family': 'serif',
    'xtick.major.size': 7,
    'xtick.major.width': 1,
    'xtick.minor.size': 4,
    'xtick.minor.width': 1,
    'ytick.major.size': 7,
    'ytick.major.width': 1,
    'ytick.minor.size': 4,
    'ytick.minor.width': 1,
    'axes.linewidth': 1,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'xtick.direction': 'in',
    'ytick.direction': 'in'
})

# print(len(err_norm))
tau_threshold = 10.
id_tau_min = np.where(tau_min>=tau_threshold)[0]
id_tau_median = np.where(tau_best>=tau_threshold)[0]
id_tau_max = np.where(tau_max>=tau_threshold)[0]
# id_x = np.where(wave[id_tau_median])[0]
print("The wave values for min :",wave_rest_fixed[id_tau_min])
print("The wave values for median :",wave_rest_fixed[id_tau_median])
print("The wave values for max :",wave_rest_fixed[id_tau_max])
# Plotting the stacked spectra
fig, axs = plt.subplots(1, 2, figsize=(15, 6))

# Left subplot (median)
axs[0].step(wave_rest_fixed, stacked_median, where='mid', color='k', lw=2.0)
axs[0].errorbar(wave_rest_fixed, stacked_median, yerr=std_median, fmt='none', ecolor='k', capsize=7)
# axs[0].errorbar(wave_rest_fixed, stacked_median, yerr=err_norm, fmt='none', ecolor='k', capsize=7)

axs[0].plot(wave_rest_fixed, flux_rest_min*6.9, color='g', linestyle='--', lw=1.8)
axs[0].plot(wave_rest_fixed, flux_rest_med*6.9, color='r', linestyle='--', lw=1.8)
axs[0].plot(wave_rest_fixed, flux_rest_max*6.9, color='g', linestyle='--', lw=1.8)
axs[0].plot(wave_rest_fixed, ermedian, 'k', linestyle='dashed', lw=2.0)
axs[0].vlines(x=1215.6701, ymin=np.min(stacked_median), ymax=0.1, colors='k', linestyles='dashed', linewidth=1.8)
axs[0].hlines(y = np.min(stacked_median)-0.05,xmin=np.min(wave_rest_fixed[id_tau_min]),
              xmax = np.max(wave_rest_fixed[id_tau_median]),color='b',linewidth=2.0,capstyle='projecting')
# axs[0].set_ylim(np.min(stacked_median)-0.05, 0.20)
# axs[0].set_ylim(-0.1, 0.20)
# axs[0].set_xlim(lo,l1)
axs[0].set_xlim(1210.,1220.)

axs[0].set_xlabel('DLA rest-frame wavelength',fontsize=20.)
axs[0].set_ylabel(r'Luminosity (10$^{40}$ erg.s$^{-1}$) ',fontsize=20)
# axs[0].set_ylim(np.min(stacked_median)-0.06, 0.25)
axs[0].set_title('Stacked Median')

# for i in range(len(wave_rest_fixed) - 1):
#     if tau[i] >= tau_threshold:
#         axs[0].fill_between([wave_rest_fixed[i], wave_rest_fixed[i+1]], -0.03, 0.1, color='orange', alpha=0.35)


# Plot blue line just below the error bar plot
# for i in range(len(wave_rest_fixed) - 1):
#     if tau[i] > tau_threshold:
        # axs[0].plot([wave_rest_fixed[i], wave_rest_fixed[i+1]], [stacked_median[i] - std_median[i], stacked_median[i+1] - std_median[i+1]], 
        #             color='blue', lw=3)

# Right subplot (mean)
axs[1].step(wave_rest_fixed, stacked_mean, where='mid', color='k', lw=2.0)
axs[1].errorbar(wave_rest_fixed, stacked_mean, yerr=std_mean, fmt='none', ecolor='k', capsize=7)
# axs[1].errorbar(wave_rest_fixed, stacked_mean, yerr=err_norm, fmt='none', ecolor='k', capsize=7)
axs[1].plot(wave_rest_fixed, flux_rest_min*6.9, color='g', linestyle='--', lw=1.8)
axs[1].plot(wave_rest_fixed, flux_rest_med*6.9, color='r', linestyle='--', lw=1.8)
axs[1].plot(wave_rest_fixed, flux_rest_max*6.9, color='g', linestyle='--', lw=1.8)
axs[1].plot(wave_rest_fixed, errmean, 'k', linestyle='dashed', lw=1.5)
axs[1].vlines(x=1215.6701, ymin=np.min(stacked_mean), ymax=0.1, colors='k', linestyles='dashed', linewidth=1.8)

axs[1].hlines(y = np.min(stacked_mean)-0.05,xmin=np.min(wave_rest_fixed[id_tau_min]),
              xmax = np.max(wave_rest_fixed[id_tau_median]),color='b',linewidth=2.0,capstyle='projecting')

axs[1].set_xlim(1210.0, 1220.0)
# axs[1].set_xlim(lo,l1)     
      
# axs[1].set_ylim(np.min(stacked_mean)-0.05, 0.20)
# axs[1].set_ylim(-0.1, 0.20)
axs[1].set_xlabel('DLA rest-frame wavelength',fontsize=20)
axs[1].set_ylabel(r'Luminosity (10$^{40}$ erg.s$^{-1}$) ',fontsize=20)
axs[1].set_ylim(np.min(stacked_mean)-0.06, 0.25)
axs[1].set_title('Weighted Mean')

plt.tight_layout()
plt.savefig('weighted_mean_stacked_figure_subplots.pdf')
# plt.savefig('significance_level_snr.pdf')
plt.show()

# ####### -------------------Full Stacked median  ----------------------###############
# print(np.array(stacked_median))
# id_w = np.where((wave_rest_fixed >= 1210.0) & (wave_rest_fixed <= 1225.0))[0]
stacked_median = np.array(stacked_median)
std_median = np.array(std_median)
ermedian = np.asarray(ermedian)
err_norm_median = np.array(err_norm)
upper_bound = np.array(upper_bound)
lower_bound = np.array(lower_bound)
best_nh1 = np.array(best_nh1)

df = pd.DataFrame({'Wavelength': wave_rest_fixed,
                   'Flux': stacked_median,
                   'Error': std_median,
                   'Error_norm':err_norm_median,
                    'Errmedian':ermedian,
                    'NHI_up': upper_bound,
                   'NHI_low':lower_bound,
                   'NHI_best':best_nh1})
###Save DataFrame to CSV
df.to_csv('stacked_median_data.csv', index=False)

# ###------------ Full Stacked Mean ------------------#######################
stacked_mean = np.array(stacked_mean)
std_mean = np.array(std_mean)
errmean = np.asarray(errmean)
err_norm_mean = np.array(err_norm)
upper_bound = np.array(upper_bound)
lower_bound = np.array(lower_bound)
best_nh1 = np.array(best_nh1)
df_mean = pd.DataFrame({'Wavelength': wave_rest_fixed,
                   'Flux': stacked_mean,
                   'Error': std_mean,
                   'Error_norm':err_norm_mean,
                   'Errmean': errmean,
                   'NHI_up': upper_bound,
                   'NHI_low':lower_bound,
                   'NHI_best':best_nh1})

###Save DataFrame to CSV
df_mean.to_csv('weighted_mean_stacked_data.csv', index=False)




# # ####### ------------------- Stacked median  ----------------------###############
# id_tau = np.where(tau_min >= tau_threshold)[0]
id_x = np.where((wave_rest_fixed > 1205.0) & (wave_rest_fixed < 1225.0))[0]
# print(wave_rest_fixed[id_tau])
# print(len(wave_rest_fixed[id_x]))
# print(len(stacked_median[id_x]))
# print(len(std_median[id_x]))
# print(len(ermedian[id_x]))
# print(len(err_norm_median[id_x]))

id_t = np.where(tau_min[id_x] <= tau_threshold)
# print(len(id_t))

######## Generate random numbers and replace values in stacked
# max_std = np.max(std_median[id_x][id_t])/2
# min_std = -0.5*(np.min(std_median[id_x][id_t])/2)
# print('the max std is',max_std,min_std)
##### Replace values in stacked with random numbers within max_std range
for idt in id_t:
    new_value = np.random.uniform(-0.0007,0.0007)
    stacked_median[id_x[idt]] = np.clip(new_value, 0.,stacked_median[id_x[idt]])  # Ensure new_value doesn't exceed original value
    # print(new_value,    stacked[id_x[idt]])  # Ensure new_value doesn't exceed original value
###### Create a DataFrame
df = pd.DataFrame({'Wavelength': wave_rest_fixed[id_x],
                   'Flux': stacked_median[id_x],
                   'Error': std_median[id_x],
                    'Errmedian':ermedian[id_x],
                   'Error_norm':err_norm_median[id_x]})
###Save DataFrame to CSV
df.to_csv('stacked_median_with_different_data.csv', index=False)


# ############------------- Stacked Weighted Mean ---------------- #################
# id_tau_mean = np.where(tau_min >= tau_threshold)[0]
id_x_mean = np.where((wave_rest_fixed > 1205.0) & (wave_rest_fixed < 1225.0))[0]
# print(wave_rest_fixed[id_x_mean])
id_t_mean = np.where(tau_min[id_x_mean] <= tau_threshold)
# print(len(id_t_mean))
##### Replace values in stacked with random numbers within max_std range
for idt_mean in id_t_mean:
    new_value_mean = np.random.uniform(-0.0007,0.0007)
    stacked_mean[id_x_mean[idt_mean]] = np.clip(new_value_mean, 0.,stacked_mean[id_x_mean[idt_mean]])  # Ensure new_value doesn't exceed original value
    

# print(stacked_mean[id_x_mean])

###### Create a DataFrame
df_mean = pd.DataFrame({'Wavelength': wave_rest_fixed[id_x_mean],        
                   'Flux': stacked_mean[id_x_mean],
                   'Error': std_mean[id_x_mean],
                   'ErrMean': errmean[id_x_mean],
                   'Error_norm':err_norm_mean[id_x_mean]})

###Save DataFrame to CSV
df_mean.to_csv('weighted_mean_stacked_different_data.csv', index=False)

# import numpy as np
# from scipy.optimize import curve_fit

# # Assuming you have the stacked spectrum stored in the variable 'stacked'
# # and the wavelength values stored in the variable 'wave_rest_fixed'

# # Define a Gaussian function
# def gaussian(x, mean, stddev):
#     amplitude = 1 / (stddev * np.sqrt(2 * np.pi))
#     return amplitude * np.exp(-0.5*(((x - mean) / stddev) ** 2) )

# # Find the maximum value of the stacked spectrum
# max_value = np.max(stacked_mean)

# # Determine the half-maximum value
# half_max = max_value / 2

# # Find the indices where the spectrum crosses the half-maximum value
# indices_above = np.where(stacked_mean >= half_max)[0]
# indices_below = np.where(stacked_mean < half_max)[0]

# # Find the wavelength values corresponding to these indices
# wavelengths_above = wave_rest_fixed[indices_above]
# wavelengths_below = wave_rest_fixed[indices_below]
# a = 1
# sigma = 2

# # Fit Gaussian curves to the spectrum
# popt_above, _ = curve_fit(gaussian, wavelengths_above, stacked_mean[indices_above], p0=[ np.mean(wavelengths_above),sigma], maxfev = 800000)
# popt_below, _ = curve_fit(gaussian, wavelengths_below, stacked_mean[indices_below], p0=[np.mean(wavelengths_below),sigma], maxfev = 800000)

# # Calculate the FWHM as the difference between the two standard deviations
# fwhm = 2.355 * np.abs(popt_above[1] - popt_below[1])

# print("Full Width at Half Maximum (FWHM):", fwhm)



# import numpy as np
# from scipy.optimize import curve_fit

# # Assuming you have the stacked spectrum stored in the variable 'stacked'
# # and the wavelength values stored in the variable 'wave_rest_fixed'

# # Define a Gaussian function
# def gaussian(x, mean, stddev):
#     amplitude = 1 / (stddev * np.sqrt(2 * np.pi))
#     return amplitude * np.exp(-0.5*(((x - mean) / stddev) ** 2) )

# # Find the maximum value of the stacked spectrum
# max_value = np.max(stacked)

# # Determine the half-maximum value
# half_max = max_value / 2

# # Find the indices where the spectrum crosses the half-maximum value
# indices_above = np.where(stacked >= half_max)[0]
# indices_below = np.where(stacked < half_max)[0]

# # Find the wavelength values corresponding to these indices
# wavelengths_above = wave_rest_fixed[indices_above]
# wavelengths_below = wave_rest_fixed[indices_below]
# a = 1
# sigma = 2

# # Fit Gaussian curves to the spectrum
# popt_above, _ = curve_fit(gaussian, wavelengths_above, stacked[indices_above], p0=[ np.mean(wavelengths_above),sigma], maxfev = 800000)
# popt_below, _ = curve_fit(gaussian, wavelengths_below, stacked[indices_below], p0=[np.mean(wavelengths_below),sigma], maxfev = 800000)

# # Calculate the FWHM as the difference between the two standard deviations
# fwhm = 2.355 * np.abs(popt_above[1] - popt_below[1])

# print("Full Width at Half Maximum (FWHM):", fwhm)


# # print(stacked)
# # print(std)

# ### Plotting start from here



# # fig = plt.figure()
# # aa = fig.add_subplot(111)
# # aa.minorticks_on()

# # # aa.fill_between(grid, 0, stacked, color='0.5', zorder=0)
# # # for i in indiv:
# # #     aa.plot(grid_indiv, i, lw=0.1, zorder=0)
# # #-0.0136
# # aa.step(wave_rest_fixed, stacked, 'k', lw=2.0, zorder=1)
# # aa.plot(wave_rest_fixed, std, 'b', lw=2.0, zorder=1)
# # aa.plot(wave_rest_fixed,errmean, 'r', lw=2.5, zorder=1)
# # aa.vlines(x=1215.67,ymin=np.min(stacked) - 0.1,ymax=np.max(stacked) + 0.1,colors='k',linestyles='dotted',linewidth=2.0)
# # aa.plot(wave_rest_fixed,flux_rest,'g',lw=1.5, zorder=1)
# # # Plotting the step plot with error bars
# # # aa.step(wave_rest_fixed, stacked, 'k', lw=0.5, zorder=1, label='Stacked Spectrum')

# # # aa.errorbar(wave_rest_fixed,stacked,yerr=std,c='k')
# # # aa.plot([wave_rest_fixed[0], wave_rest_fixed[-1]], [0, 0], lw=1, ls='--', zorder=1, color='k')
# # # aa.plot(wave_rest_fixed, std, 'b', lw=0.5, zorder=1, label=f'$\sigma$')
# # # aa.plot(wave_rest_fixed,errmean, 'r', lw=0.5, zorder=1, label='Error on the mean')
# # # aa.vlines(x=1215.67,ymin=-0.1,ymax=0.2,colors='k',linestyles='dotted')

# # # aa.plot(wave_rest_fixed,flux_rest,'g',lw=0.5, zorder=1, label=f'Voigt Profile_{value_1}')
# # for i in range(len(wave_rest_fixed) - 1):
# #     if tau[i] > tau_threshold:
# #         aa.fill_between([wave_rest_fixed[i], wave_rest_fixed[i+1]], -0.1, 0.1, color='orange', alpha=0.35)
# #         # aa.hlines(y = np.min(stacked)-0.05,xmin=wave_rest_fixed[i], xmax= wave_rest_fixed[i+1],lw = 2.0, color='blue', alpha=0.5)


# # aa.legend()
# # aa.set_xlim(wave_rest_fixed[0], wave_rest_fixed[-1])
# # aa.set_xlim(1200.0, 1225.)
# aa.set_ylim(np.min(stacked)-0.11,np.max(stacked)+0.11)

# aa.set_title('Median',fontsize=mp.rcParams['xtick.labelsize'])
# # aa.set_xlabel('$\lambda$ ($\mathrm{\AA}$)',fontsize=mp.rcParams['ytick.labelsize'])
# # aa.set_ylabel('Normalised Flux Density [unitless]')
# aa.set_xlabel("DLA's rest-frame wavelength ($\mathrm{\AA}$)", fontsize=mp.rcParams['xtick.labelsize'])
# aa.set_ylabel('Normalised Flux', fontsize=mp.rcParams['ytick.labelsize'])
# plt.savefig('median_stacked_figure.pdf')
# plt.show()

# # # Define conditions for replacement
# # id_tau = np.where(tau >= tau_threshold)[0]
# # id_x = np.where((wave_rest_fixed > 1210.0) & (wave_rest_fixed < 1220.0))[0]
# # print(wave_rest_fixed[id_x])

# # id_t = np.where(tau[id_x] < 10)
# # print(len(id_t))
# # # Generate random numbers and replace values in stacked
# # max_std = np.max(std[id_x][id_t])/2
# # min_std = -1*(np.min(std[id_x][id_t])/2)
# # print('the max std is',max_std,min_std)
# # # Replace values in stacked with random numbers within max_std range
# # for idt in id_t[0]:
# #     new_value = np.random.uniform(min_std, max_std)
# #     stacked[id_x[idt]] = np.clip(new_value, 0, stacked[id_x[idt]])  # Ensure new_value doesn't exceed original value
# #     # print(new_value,    stacked[id_x[idt]])  # Ensure new_value doesn't exceed original value
# # # # Create a DataFrame
# # df = pd.DataFrame({'Wavelength': wave_rest_fixed[id_x],
# #                    'Flux': stacked[id_x],
# #                    'Error': std[id_x]})

# # Save DataFrame to CSV
# # df.to_csv('integrated_data.csv', index=False)

# # plot(wave_rest_fixed, stacked,std, errmean,indiv=rebinned,grid_indiv=wave_rest_fixed,nh1 = value_1,cnr =value_2,sigma=7)

# # rebinned = specstack.regrid(file,wave_rest_fixed)
