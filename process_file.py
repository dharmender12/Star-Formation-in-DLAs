# # import os
# # from astropy.io import fits
# # import VoigtFit
# # # from line_components import profile_components,modify_components
# # from line_choice import velocity_plot


# # def process_file(filename,id, wave, norm_flux,z_cnn, res_UVES, norm_err, logNHI, increment_b=0, increment_N=0, modify_params=True,max_attempts=15):
    
# #     line_choices, centroid_value, primary_lines, primary_components, primary_line_values,z_sys = velocity_plot(filename, id, wave, norm_flux, z_cnn)

# #     # initial_components,z_sys = profile_components(wave,norm_flux,z_cnn,centroid_value, primary_lines, primary_components, primary_line_values,z_sys, increment_b=0, increment_N=0)
# #     # print('the initial components are :',initial_components[0][0])
    
# #     z_correct = int(input("Is the redshift correct or not (choose 1, 2, or 3): "))
    
# #     dataset = VoigtFit.DataSet(z_sys)
# #     dataset.add_data(wave, norm_flux, 299792. / res_UVES, err=norm_err, normalized=True)

# #     # Add lines
# #     dataset.add_line(primary_lines[0], velspan=(-1500, 1500))
# #     for line in primary_lines:
# #         if line != primary_lines[0]:
# #             dataset.add_line(line, velspan=(-1500, 1500))


# #     # Initial guesses
# #     b_values = [26.0, 29.0,38.0,31.0]
# #     logN_values = [12.5, 14.8,17.5,14.2]
        
# #     def add_components(dataset, z_sys, primary_components, b_values, logN_values):
# #         dataset.add_component(primary_components[0], z_sys, b_values[1], logN_values[1], var_z=1)
# #         dataset.add_component(primary_components[0], z_sys, b_values[2], logN_values[2], var_z=1)
# #         dataset.add_component(primary_components[1], z_sys, b_values[0], logN_values[0], tie_z=f'z0_{primary_components[0]}')
# #         dataset.add_component(primary_components[1], z_sys, b_values[3] , logN_values[3], tie_z=f'z1_{primary_components[0]}')

    
# #     add_components(dataset, z_sys, primary_components, b_values, logN_values)

# #     def handle_fitting_and_plotting():
# #         attempt = 0
# #         while attempt < max_attempts:
# #             try:
# #                 dataset.prepare_dataset(norm=False)
# #                 dataset.prepare_dataset()
# #                 popt, chi2 = dataset.fit()
# #                 print("The chi2 is :",chi2)
# #                 dataset.plot_fit(filename=f'./vp_fit/{os.path.basename(filename)}')

# #                 if logNHI:
# #                     dataset.print_metallicity(*logNHI)

# #                 dataset.print_total()
# #                 dataset.save_parameters(f'./vp_fit/{os.path.basename(filename)}.pars')
# #                 return True

# #             except Exception as e:
# #                 print(f"Attempt {attempt + 1} failed: {e}")
# #                 attempt += 1

# #         print("Failed to fit after maximum attempts.")
# #         return False

# #     if not handle_fitting_and_plotting() and modify_params:
# #         # Clear previous components and try modifying parameters
# #         # dataset.clear_components()
# #         # modified_components,z_sys = modify_components(wave,norm_flux,centroid_value, z_sys,increment_b=0.5, increment_N=0.6)
# #         # modified_components = modify_components(primary_line, primary_component, z_cnn, increment_b, increment_N)
# #         b_incremented = [28.0 + increment_b, 31.0 + increment_b,38.0 - increment_b ,24.0 + increment_b]
# #         logN_incremented = [13.5 + increment_N, 14.8 + increment_N,17.5 - increment_N,13.2 + increment_N]
# #         add_components(dataset, z_sys, primary_components, b_incremented, logN_incremented)

# #         if not handle_fitting_and_plotting():
# #             print("Failed to fit even after modifying parameters.")

# #     return dataset,z_correct,z_sys,primary_lines,centroid_value


# import os
# from astropy.io import fits
# import VoigtFit
# # from line_components import profile_components, modify_components
# from line_choice import velocity_plot


# def process_file(filename, id, wave, norm_flux, z_cnn, res_UVES, norm_err, logNHI, increment_b=0, increment_N=0, modify_params=True, max_attempts=15):
    
#     line_choices, centroid_value, primary_lines, primary_components, primary_line_values, z_sys = velocity_plot(filename, id, wave, norm_flux, z_cnn)

#     # initial_components, z_sys = profile_components(wave, norm_flux, z_cnn, centroid_value, primary_lines, primary_components, primary_line_values, z_sys, increment_b=0, increment_N=0)
#     # print('the initial components are :', initial_components[0][0])
    
#     z_correct = int(input("Is the redshift correct or not (choose 1, 2, or 3): "))
    
#     dataset = VoigtFit.DataSet(z_sys)
#     dataset.add_data(wave, norm_flux, 299792. / res_UVES, err=norm_err, normalized=True)

#     # Add lines
#     dataset.add_line(primary_lines[0], velspan=(-1500, 1500))
#     for line in primary_lines:
#         if line != primary_lines[0]:
#             dataset.add_line(line, velspan=(-1500, 1500))

#     # Initial guesses
#     b_values = [26.0, 29.0, 38.0, 31.0]
#     logN_values = [12.5, 14.8, 17.5, 14.2]
        
#     def add_components(dataset, z_sys, primary_components, b_values, logN_values):
#         dataset.add_component(primary_components[0], z_sys, b_values[1], logN_values[1], var_z=1)
#         dataset.add_component(primary_components[0], z_sys, b_values[2], logN_values[2], var_z=1)
#         dataset.add_component(primary_components[1], z_sys, b_values[0], logN_values[0], tie_z=f'z0_{primary_components[0]}')
#         dataset.add_component(primary_components[1], z_sys, b_values[3], logN_values[3], tie_z=f'z1_{primary_components[0]}')

#     add_components(dataset, z_sys, primary_components, b_values, logN_values)

#     def handle_fitting_and_plotting():
#         attempt = 0
#         last_chi2 = None
#         same_chi2_count = 0
#         while attempt < max_attempts:
#             try:
#                 dataset.prepare_dataset(norm=False)
#                 dataset.prepare_dataset()
#                 popt, chi2 = dataset.fit()
#                 print("The chi2 is:", chi2)
                
#                 if last_chi2 is not None and abs(last_chi2 - chi2) < 1e-6:
#                     same_chi2_count += 1
#                 else:
#                     same_chi2_count = 0

#                 if same_chi2_count > 1:  # If chi2 is the same for more than one attempt
#                     print("Chi2 did not improve. Saving data and moving to the next file.")
#                     # dataset.plot_fit(filename=f'./vp_fit/{os.path.basename(filename)}')
#                     break

#                 last_chi2 = chi2
                
#                 dataset.plot_fit(filename=f'./vp_fit/{os.path.basename(filename)}')

#                 if logNHI:
#                     dataset.print_metallicity(*logNHI)

#                 dataset.print_total()
#                 dataset.save_parameters(f'./vp_fit/{os.path.basename(filename)}.pars')
#                 return True

#             except Exception as e:
#                 print(f"Attempt {attempt + 1} failed: {e}")
#                 attempt += 1

#         print("Failed to fit after maximum attempts.")
#         return False

#     if not handle_fitting_and_plotting() and modify_params:
#         # Clear previous components and try modifying parameters
#         # dataset.clear_components()
#         # modified_components, z_sys = modify_components(wave, norm_flux, centroid_value, z_sys, increment_b=0.5, increment_N=0.6)
#         # modified_components = modify_components(primary_line, primary_component, z_cnn, increment_b, increment_N)
#         b_incremented = [26.0 + increment_b, 32.0 + increment_b, 39.0 - increment_b, 24.0 + increment_b]
#         logN_incremented = [13.5 + increment_N, 14.8 + increment_N, 17.5 - increment_N, 13.2 + increment_N]
#         add_components(dataset, z_sys, primary_components, b_incremented, logN_incremented)

#         if not handle_fitting_and_plotting():
#             print("Failed to fit even after modifying parameters.")

#     return dataset, z_correct, z_sys, primary_lines, centroid_value


import os
from astropy.io import fits
import VoigtFit
# from line_components import profile_components, modify_components
from line_choice import velocity_plot
import numpy as np
from scipy.optimize import minimize



def loss_function(params, x, y):
    mu, sigma, lambda_0, z = params
    lambda_obs = lambda_0 * (1 + z)
    amp = 1 / np.sqrt(2 * np.pi) * sigma
    model = -amp * np.exp(-((x - lambda_obs) ** 2) / (sigma ** 2))
    loss = np.sum((y - model) ** 2)  # Sum of squared residuals
    return loss

def process_file(filename, id, wave, norm_flux, z_cnn, res_UVES, norm_err, logNHI, increment_b=0, increment_N=0, modify_params=True, max_attempts=15):
    
    line_choices, centroid_value, primary_lines, primary_components, primary_line_values, z_sys = velocity_plot(filename, id, wave, norm_flux, z_cnn)

    # initial_components, z_sys = profile_components(wave, norm_flux, z_cnn, centroid_value, primary_lines, primary_components, primary_line_values, z_sys, increment_b=0, increment_N=0)
    # print('the initial components are :', initial_components[0][0])
    
    z_correct = int(input("Is the redshift correct or not (choose 1, 2, or 3): "))

    if z_correct == 1:
        z_used = z_cnn
    else:
        z_used = z_sys
    
    dataset = VoigtFit.DataSet(z_used)
    dataset.add_data(wave, norm_flux, 299792. / res_UVES, err=norm_err, normalized=True)

    # Add lines
    dataset.add_line(primary_lines[0], velspan=(-1500, 1500))
    for line in primary_lines:
        if line != primary_lines[0]:
            dataset.add_line(line, velspan=(-1500, 1500))

    # Initial guesses
    b_values = [26.0, 29.0, 38.0, 31.0]
    logN_values = [12.5, 14.8, 17.5, 14.2]
    # Define the loss function (to be minimized)
    initial_guess = [26.9, 14.3, centroid_value, z_sys]

    # Bounds for b and logN
    bounds = [(10.0, 40.0), (11.5, 18.9), (centroid_value - 1.0, centroid_value + 1.0), (z_sys - 0.0002, z_sys + 0.0002)]
    
    # Perform minimization
    result = minimize(loss_function, initial_guess, args=(wave, norm_flux), bounds=bounds)
    b_opt, logN, lambda_0, z_opt = result.x
    
    def add_components(dataset, z_used, primary_components, b_values, logN_values):
        dataset.add_component(primary_components[0], z_used, b_opt, logN, var_z=1)
        dataset.add_component(primary_components[0], z_used, b_values[1], logN_values[1], var_z=1)
        dataset.add_component(primary_components[1], z_used, b_values[0], logN_values[0], tie_z=f'z0_{primary_components[0]}')
        dataset.add_component(primary_components[1], z_used, b_values[3], logN_values[3], tie_z=f'z1_{primary_components[0]}')

    add_components(dataset, z_used, primary_components, b_values, logN_values)

    def handle_fitting_and_plotting():
        attempt = 0
        last_chi2 = None
        same_chi2_count = 0
        while attempt < max_attempts:
            try:
                dataset.prepare_dataset(norm=False)
                dataset.prepare_dataset()
                popt, chi2 = dataset.fit()
                print("The chi2 is:", chi2)
                
                if last_chi2 is not None and abs(last_chi2 - chi2) < 1e-6:
                    same_chi2_count += 1
                else:
                    same_chi2_count = 0

                if same_chi2_count > 1:  # If chi2 is the same for more than one attempt
                    print("Chi2 did not improve. Saving data and moving to the next file.")
                    # dataset.plot_fit(filename=f'./vp_fit/{os.path.basename(filename)}')
                    break

                last_chi2 = chi2
                
                dataset.plot_fit(filename=f'./vp_fit/{os.path.basename(filename)}')

                if logNHI:
                    dataset.print_metallicity(*logNHI)

                dataset.print_total()
                dataset.save_parameters(f'./vp_fit/{os.path.basename(filename)}.pars')
                return True

            except Exception as e:
                print(f"Attempt {attempt + 1} failed: {e}")
                attempt += 1

        print("Failed to fit after maximum attempts.")
        return False

    if not handle_fitting_and_plotting() and modify_params:
        # Clear previous components and try modifying parameters
        # dataset.clear_components()
        # modified_components, z_used = modify_components(wave, norm_flux, centroid_value, z_used, increment_b=0.5, increment_N=0.6)
        # modified_components = modify_components(primary_line, primary_component, z_cnn, increment_b, increment_N)
        b_incremented = [20.0 + increment_b, 26.0 + increment_b, 33.0 - increment_b, 25.0 + increment_b]
        logN_incremented = [11.6 + increment_N, 14.5 - increment_N, 17.5 - increment_N, 13.2 + increment_N]
        add_components(dataset, z_used, primary_components, b_incremented, logN_incremented)

        if not handle_fitting_and_plotting():
            print("Failed to fit even after modifying parameters.")

    return dataset, z_correct, z_sys, primary_lines, centroid_value
