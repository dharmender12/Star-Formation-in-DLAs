import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import numpy as np
import os
from matplotlib import rcParams
from astropy.io import fits
# Set the default font family
rcParams['font.family'] = 'DejaVu Sans'



# Define the inverted Gaussian function for fitting
def inverted_gaussian(x, amp, mu, sigma):
    # Prevent overflow by clipping sigma to a reasonable range
    sigma = np.clip(sigma, 1e-3, 10)
    return -amp * np.exp((x - mu)**2 / (2 * sigma**2))


def choose_lines():
    print("Choose the primary lines for the process (comma-separated):")
    print("0: HI_1215")
    print("1: CII (1334.5)")
    print("2: Al II (1670.0)")
    print("3: Si II (1526.7)")
    print("4: O I (1302.2)")
    print("5: FeII_2382")
    print("6: FeII_2344")
    print("7: FeII_2600")
    # line_choice = input(int("Enter the number for plotting"))
    choices = input("Enter the numbers fitting  (e.g., 1,3,5,6,7): ").split(',') 
    # choices = input("Enter the numbers fitting  (e.g., 1,3): ").replace(',','.').split()
    # choices = input("Enter the numbers fitting  (e.g., 1,3): ").split('.')
    line_choices = [int(choice) for choice in choices]
    
    line_dict = {
        0: 'HI_1215',
        1: 'CII_1334',
        2: 'AlII_1670',
        3: 'SiII_1526',
        4: 'OI_1302',
        5:'FeII_2382',
        6: 'FeII_2344',
        7: 'FeII_2600'
        
    }
    
    primary_line_value_dict = {
        0: 1215.6701,
        1: 1334.5,
        2: 1670.0,
        3: 1526.7,
        4: 1302.2,
        5:2382.764,
        6:2344.212,
        7:2600.172
    }
    
    component_dict = {
        0: 'HI',
        1: 'CII',
        2: 'AlII',
        3: 'SiII',
        4: 'OI',
        5: 'FeII_1',
        6: 'FeII_2',
        7: 'FeII_3'
    }
    color_dict = {
        0: 'b',
        1: 'g',
        2: 'm',
        3: 'y',
        4: 'c',
        5: 'r',
        6: 'grey',
        7: 'cyan'
    }

    primary_lines = []
    primary_components = []
    primary_line_values = []
    primary_colors = []

    for line_choice in line_choices:
        if line_choice not in line_dict:
            print(f"Invalid choice: {line_choice}. Please choose a number between 0 and 6.")
            return choose_lines()
        
        primary_lines.append(line_dict[line_choice])
        primary_components.append(component_dict[line_choice])
        primary_line_values.append(primary_line_value_dict[line_choice])
        primary_colors.append(color_dict[line_choice])
    
    
    
    return line_choices, primary_lines, primary_components, primary_line_values,primary_colors

def velocity_plot(filename, id, wave, normalised_flux, z_cnn):
    # Define wavelength lines
    dla_wave = 1215.6701
    c2_wave = 1334.5
    alii_1670 = 1670.0
    si2_wave = 1526.7
    abs_o1_1302 = 1302.2
    feii_2382 = 2382.764
    feii_2344 = 2344.21
    feii_2600 = 2600.172
    z_dla = (wave / dla_wave) - 1
    z_si2 = (wave / si2_wave) - 1
    z_cii = (wave / c2_wave) - 1
    z_o1 = (wave / abs_o1_1302) - 1
    z_alii = (wave / alii_1670) - 1
    z_feii_1 = (wave / feii_2382) - 1
    z_feii_2 = (wave / feii_2344) - 1
    z_feii_3 = (wave / feii_2600) - 1
    
    speed_of_light = 3e5  # km/s

    velocity_o1 = speed_of_light * (((1 + z_o1)**2 - (1 + z_cnn)**2) / ((1 + z_o1)**2 + (1 + z_cnn)**2))
    velocity_dla = speed_of_light * (((1 + z_dla)**2 - (1 + z_cnn)**2) / ((1 + z_dla)**2 + (1 + z_cnn)**2))
    velocity_si2 = speed_of_light * (((1 + z_si2)**2 - (1 + z_cnn)**2) / ((1 + z_si2)**2 + (1 + z_cnn)**2))
    velocity_cii = speed_of_light * (((1 + z_cii)**2 - (1 + z_cnn)**2) / ((1 + z_cii)**2 + (1 + z_cnn)**2))
    velocity_alii = speed_of_light * (((1 + z_alii)**2 - (1 + z_cnn)**2) / ((1 + z_alii)**2 + (1 + z_cnn)**2))
    velocity_feii_1 = speed_of_light * (((1 + z_feii_1)**2 - (1 + z_cnn)**2) / ((1 + z_feii_1)**2 + (1 + z_cnn)**2))
    velocity_feii_2 = speed_of_light * (((1 + z_feii_2)**2 - (1 + z_cnn)**2) / ((1 + z_feii_2)**2 + (1 + z_cnn)**2))
    velocity_feii_3 = speed_of_light * (((1 + z_feii_3)**2 - (1 + z_cnn)**2) / ((1 + z_feii_3)**2 + (1 + z_cnn)**2))
    # velocity_feii_1 = speed_of_light * (((1 + z_alii)**2 - (1 + z_cnn)**2) / ((1 + z_alii)**2 + (1 + z_cnn)**2))
    plt.figure(figsize=(12, 8))
    # Plot the velocity space
    # plt.step(velocity_dla, normalised_flux, label='HI (1215.6701)', color='b')
    plt.plot(velocity_cii, normalised_flux, label='CII (1334.5)', color='g')
    plt.plot(velocity_si2, normalised_flux, label='SiII (1526.7)', color='y')
    plt.plot(velocity_o1, normalised_flux, label='OI (1302.2)', color='c')
    plt.plot(velocity_alii, normalised_flux, label='AlII (1670.0)', color='m')
    plt.plot(velocity_feii_1, normalised_flux, label='FeII_2382', color='r')
    plt.plot(velocity_feii_2, normalised_flux, label='Fe_2344', color='orange')
    plt.plot(velocity_feii_3, normalised_flux, label='Fe_2600', color='darkgreen')
    plt.title(f'Vel plot for ID {id} at redshift z={z_cnn}')
    plt.axvline(x=0.0, color='k', linewidth=0.25)
    plt.xlim(-2000, 2000)
    plt.ylim(-0.50,2.0)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{os.path.basename(filename)}_vel.pdf')
    plt.show()
    plt.close()

    # Get user choice for spectrum plot
    line_choices, primary_lines, primary_components, primary_line_values,line_colors = choose_lines()
    
    for i in range(len(primary_line_values)): 
        # Define the region around the line for fitting
        line_center = primary_line_values[i] * (1 + z_cnn)
        fit_region = (wave > line_center - 10) & (wave < line_center + 10)
        wave_fit = wave[fit_region]
        flux_fit = normalised_flux[fit_region]

        # Improved initial guess for the Gaussian parameters
    
        centroid_value = wave_fit[np.argmin(flux_fit)]

        z_sys =( centroid_value / primary_line_values[i]) -1

        plt.figure(figsize=(12, 8))
        plt.plot(wave, normalised_flux, label="Spectrum")
        plt.axvline(x=line_center, label=f"{primary_lines[i]}_centroid_{centroid_value}", color=f'{line_colors[i]}')
        # plt.plot(wave_fit, fitted_curve, label='Inverted Gaussian Fit')
        plt.xlim(line_center - 50, line_center + 50)
        plt.ylim(-1,2.0)
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Normalized Flux')
        plt.title(f'Spectrum for ID {id} at redshift z={z_cnn}')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./spec_plots/{os.path.basename(filename)}_spec.pdf')
        plt.show()
        plt.close()
        
    return line_choices, centroid_value, primary_lines, primary_components, primary_line_values,z_sys



filename = 'dla_500kpc_237.fits'#dla_150kpc_23.fits'#'dla_100kpc_15.fits'
spectrum = fits.open(f'{filename}')
data = spectrum[1].data
wave = data['WAVE    ']
flux = data['STACK_MEDIAN ']
id = 0
z_cnn = 0

line_choices, centroid_value, primary_lines, primary_components, primary_line_values,z_sys = velocity_plot(filename, id, wave,flux, z_cnn)

# import matplotlib.pyplot as plt
# from scipy.optimize import minimize,differential_evolution,curve_fit
# import numpy as np
# import os

# # Define the inverted Gaussian function for fitting
# def inverted_gaussian(x,amp, mu, sigma):
#     # amp = 1/ np.sqrt(2 * np.pi) * sigma
#     return -amp * np.exp((x - mu)**2 / (2 * sigma**2))

# def choose_line():
#     print("Choose the primary line for the process:")
#     print("0: HI_1215")
#     print("1: C II (1334.5)")
#     print("2: Al II (1670.0)")
#     print("3: Si II (1526.7)")
#     print("4: O I (1302.2)")

#     line_choice = int(input("Enter the number corresponding to your choice: "))
    
#     line_dict = {
#         0: 'HI_1215',
#         1: 'CII_1334',
#         2: 'AlII_1670',
#         3: 'SiII_1526',
#         4: 'OI_1302.2'
#     }
    
#     primary_line_value_dict = {
#         0: 1215.6701,
#         1: 1334.5,
#         2: 1670.0,
#         3: 1526.7,
#         4: 1302.2
#     }
    
#     component_dict = {
#         0: 'HI',
#         1: 'CII',
#         2: 'AlII',
#         3: 'SiII',
#         4: 'OI'
#     }
#     color_dict = {
#         0:'b',
#         1:'g',
#         2:'y',
#         3:'r',
#         4:'m'
#     }
    
#     if line_choice not in line_dict:
#         print("Invalid choice. Please choose a number between 0 and 4.")
#         return choose_line()
    
#     primary_line = line_dict[line_choice]
#     primary_component = component_dict[line_choice]
#     primary_line_value = primary_line_value_dict[line_choice]
#     primary_color = color_dict[line_choice]
#     return line_choice, primary_line, primary_component, primary_line_value,primary_color

# def velocity_plot(filename, id, wave, normalised_flux, z_cnn):
#     # Define wavelength lines
#     dla_wave = 1215.6701
#     c2_wave = 1334.5
#     alii_1670 = 1670.0
#     si2_wave = 1526.7
#     abs_o1_1302 = 1302.2
#     z_dla = (wave / dla_wave) - 1
#     z_si2 = (wave / si2_wave) - 1
#     z_cii = (wave / c2_wave) - 1
#     z_o1 = (wave / abs_o1_1302) - 1
#     z_alii = (wave / alii_1670) - 1
#     speed_of_light = 3e5  # km/s

#     velocity_o1 = speed_of_light * (((1 + z_o1)**2 - (1 + z_cnn)**2) / ((1 + z_o1)**2 + (1 + z_cnn)**2))
#     velocity_dla = speed_of_light * (((1 + z_dla)**2 - (1 + z_cnn)**2) / ((1 + z_dla)**2 + (1 + z_cnn)**2))
#     velocity_si2 = speed_of_light * (((1 + z_si2)**2 - (1 + z_cnn)**2) / ((1 + z_si2)**2 + (1 + z_cnn)**2))
#     velocity_cii = speed_of_light * (((1 + z_cii)**2 - (1 + z_cnn)**2) / ((1 + z_cii)**2 + (1 + z_cnn)**2))
#     velocity_alii = speed_of_light * (((1 + z_alii)**2 - (1 + z_cnn)**2) / ((1 + z_alii)**2 + (1 + z_cnn)**2))

#     plt.figure(figsize=(12, 8))
#     # Plot the velocity space
#     plt.step(velocity_dla, normalised_flux, label='HI (1215.6701)', color='b')
#     plt.plot(velocity_cii, normalised_flux, label='CII (1334.5)', color='g')
#     plt.plot(velocity_si2, normalised_flux, label='SiII (1526.7)', color='y')
#     plt.plot(velocity_o1, normalised_flux, label='OI (1302.2)', color='c')
#     plt.plot(velocity_alii, normalised_flux, label='AlII (1670.0)', color='m')
#     plt.title(f'Vel plot for ID {id} at redshift z={z_cnn}')
#     plt.axvline(x=0.0, color='k', linewidth=0.25)
#     plt.xlim(-2000, 2000)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(f'./vel_plots/{os.path.basename(filename)}_vel.pdf')
#     plt.show()
#     plt.close()

#     # Get user choice for spectrum plot
#     line_choice, primary_line, primary_component, primary_line_value,primary_color = choose_line()
    
#     # Define the region around the line for fitting
#     line_center = primary_line_value * (1 + z_cnn)
#     fit_region = (wave > line_center - 10) & (wave < line_center + 10)
#     wave_fit = wave[fit_region]
#     flux_fit = normalised_flux[fit_region]

#     # Initial guess for the Gaussian parameters
#     amp_guess = np.max(flux_fit) - np.min(flux_fit)
#     mu_guess = np.mean(wave_fit)
#     sigma_guess = np.std(wave_fit)
#     initial_guess = [amp_guess,mu_guess, sigma_guess]
    
#     def objective(params):
#         return np.sum((inverted_gaussian(wave_fit, *params) - flux_fit)**2)
    
#     bounds = [(0.01, 1.5), (np.min(wave_fit)+1, np.max(wave_fit)-1), (0.001, 6.0)]

#     try:
#         # Fit the Gaussian using minimize
#         result = minimize(objective, initial_guess, bounds=bounds)
#         # Fit the Gaussian
#         # popt, _ = curve_fit(inverted_gaussian, wave_fit, flux_fit, p0=initial_guess)#, bounds=bounds)
#         popt = result.x
#         centroid_value = popt[1]
#         print("the centroid values are: ",centroid_value)
#         fitted_curve = inverted_gaussian(wave_fit, *popt)
#     except RuntimeError:
#         print(f"Fit failed for {filename}, using line_center as centroid")
#         centroid_value = line_center
#         popt = [line_center, sigma_guess]
#         fitted_curve = inverted_gaussian(wave_fit, *popt)
    
#     plt.figure(figsize=(12, 8))
#     plt.plot(wave, normalised_flux, label="Spectrum")
#     plt.axvline(x=(1 + z_cnn) * primary_line_value, label=f"{primary_line}", color=f'{primary_color}')
#     plt.plot(wave_fit, fitted_curve, label='Inverted Gaussian Fit')
#     plt.xlim((1 + z_cnn) * primary_line_value - 50, (1 + z_cnn) * primary_line_value + 50)
#     plt.xlabel('Wavelength (Å)')
#     plt.ylabel('Normalized Flux')
#     plt.title(f'Spectrum for ID {id} at redshift z={z_cnn} and centroid_value {centroid_value}')
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(f'./spec_plots/{os.path.basename(filename)}_spec.pdf')
#     plt.show()
#     plt.close()
    
#     return line_choice, centroid_value, primary_line, primary_component, primary_line_value


