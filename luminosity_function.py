import Lya_zelda as Lya
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18
import astropy.units as u
from scipy.optimize import curve_fit

z_avg = 2.679

def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

def compute_snr(wavelength, flux, error, wave_range):
    # Extract relevant region based on wavelength range
    region = np.where((wavelength >= wave_range[0]) & (wavelength <= wave_range[1]))[0]
    flux_region = flux[region]
    error_region = error[region]

    # Integrated flux (signal)
    integrated_flux = np.trapz(flux_region, wavelength[region])

    # Noise estimation (square root of the sum of squared errors)
    noise_level = np.sqrt(np.trapz(error_region,wavelength[region])**2)

    # SNR calculation
    snr = integrated_flux / noise_level
    
    return snr

def sfr(data):

    wave = np.array(data['Wavelength'])
    flux = np.array(data['Flux'])
    err = np.array(data['Error_norm'])
    id_x = np.where((wave >= 1213.4833) & (wave <= 1217.8027))[0]
    id_tau = np.where((wave >= 1215.6700) & (wave <= 1217.8027))[0]

    print("wave:", wave)
    w_Arr = wave[id_x]
    f_Arr = flux[id_x]
    s_Arr = err[id_x]

    integrated_flux = np.trapz(f_Arr, w_Arr)

    w_Arr_tau = wave[id_tau]
    f_Arr_tau = flux[id_tau]
    s_Arr_tau = err[id_tau]

    w_Arr_tau_cm = w_Arr_tau
    integrated_flux_tau_cm = np.trapz(f_Arr_tau, w_Arr_tau_cm)
    integrated_flux_tau_cm_all = np.trapz(f_Arr, w_Arr)

    integrated_err = np.trapz(s_Arr_tau,w_Arr_tau)
    integrated_err_all = np.trapz(s_Arr,w_Arr)

    lum_dist = Planck18.luminosity_distance(z_avg)
    lum_dist_cm = lum_dist.to(u.cm).value

    luminosity = 4 * np.pi * (lum_dist_cm ** 2) * integrated_flux_tau_cm * (1e-17) # erg/s
    luminosity_all = 4 * np.pi * (lum_dist_cm ** 2) * integrated_flux_tau_cm_all * (1e-17)

    sfr = 9.08 * 10**(-43) * luminosity

    known_wavelength = 1215.67
    speed_of_light = 3e5 # km/s

    popt, pcov = curve_fit(gaussian, w_Arr, f_Arr, p0=[max(flux), wave[np.argmax(flux)], 1])

    A, mu, sigma = popt
    FWHM = 2 * np.sqrt(2 * np.log(2)) * sigma
    FWHM_km_s = (FWHM / known_wavelength) * speed_of_light  # Convert FWHM to km/s
    velocity_offset = (mu - known_wavelength) / known_wavelength * speed_of_light

    # Calculate the uncertainties for the parameters
    perr = np.sqrt(np.diag(pcov))
    A_err, mu_err, sigma_err = perr

    # Propagate the error for the velocity offset
    velocity_offset_err = np.sqrt(
        (speed_of_light / known_wavelength) ** 2 * (mu_err ** 2)
    )
    
    plt.plot(wave, flux, label='Data 1')
    plt.plot(wave, gaussian(wave, *popt), label='Gaussian Fit 1')
    plt.legend()
    plt.show()
    plt.close()

    popt, pcov = curve_fit(gaussian,w_Arr_tau , f_Arr_tau, p0=[max(flux), wave[np.argmax(flux)], 1])

    A_1, mu_1, sigma_1 = popt

    # Use the new SNR formula
    wave_range = (1215.6701, 1217.8027)
    SNR = compute_snr(wave, flux, err, wave_range)
    snr_all = compute_snr(wave, flux, err, (1213.4835, 1217.8027))

    print(f"Data - Sigma: {sigma}, SNR: {SNR}, SNR All: {snr_all}, FWHM: {FWHM} Å, FWHM (km/s): {FWHM_km_s}, Velocity Offset: {velocity_offset} ± {velocity_offset_err}")
    print(f'Luminosity: {luminosity}, SFR: {sfr}, and Luminosity Full: {luminosity_all}')
    return sfr, luminosity, luminosity_all

data_1 = pd.read_csv('stacked_median_with_different_data.csv')
data_2 = pd.read_csv('weighted_mean_stacked_different_data.csv')

sfr_med, lum_med, lum_med_all = sfr(data_1)
with open('median_sfr_data.txt', 'w') as f:
    f.write(f'sfr_med: {sfr_med}\n')
    f.write(f'lum_med: {lum_med}\n')
    f.write(f'lum_med_all: {lum_med_all}\n')

sfr_mean, lum_mean, lum_mean_all = sfr(data_2)
with open('mean_sfr_data.txt', 'w') as f:
    f.write(f'sfr_mean: {sfr_mean}\n')
    f.write(f'luminosity_mean: {lum_mean}\n')
    f.write(f'luminosity_mean_all: {lum_mean_all}\n')
