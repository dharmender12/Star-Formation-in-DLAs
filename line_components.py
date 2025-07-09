import numpy as np
from scipy.optimize import minimize, curve_fit
from line_choice import velocity_plot

# Define the Gaussian profile function
def gaussian_profile(x, amp, mu, sigma, lambda_0, z):
    lambda_obs = lambda_0 * (1 + z)
    return amp * np.exp(-((x - lambda_obs) ** 2) / (2 * sigma ** 2))

# Define the loss function (to be minimized)
def loss_function(params, x, y):
    mu, sigma, lambda_0, z = params
    lambda_obs = lambda_0 * (1 + z)
    amp = 1 / np.sqrt(2 * np.pi) * sigma
    model = -amp * np.exp(-((x - lambda_obs) ** 2) / (sigma ** 2))
    loss = np.sum((y - model) ** 2)  # Sum of squared residuals
    return loss

# Initial components
def profile_components(wave, norm_flux, centroid_value, z_sys,increment_b=0,increment_N=0):
    initial_components = []
    # Initial guess for parameters (b, logN, lambda_0, z)
    initial_guess = [16.9, 14.3, centroid_value, z_sys]

    # Bounds for b and logN
    bounds = [(10.0, 40.0), (11.5, 18.9), (centroid_value - 1.0, centroid_value + 1.0), (z_sys - 0.0002, z_sys + 0.0002)]
    
    # Perform minimization
    result = minimize(loss_function, initial_guess, args=(wave, norm_flux), bounds=bounds)
    
    # Extract optimized parameters
    b_opt, logN, lambda_0, z_opt = result.x
    initial_components = [b_opt, logN, lambda_0, z_opt]
    print("optimized b: ", b_opt)
    
    return initial_components, z_sys

def modify_components(wave, norm_flux, centroid_value, z_sys, increment_b=0.5, increment_N=0.6):
    modified_components = []

    initial_guess = [25.0, 15.4, centroid_value, z_sys]
    # Define bounds for the parameters
    bounds = ([1.0, 11.5, centroid_value - 1.0, z_sys - 0.0002],
              [10.0, 18.9, centroid_value + 1.0, z_sys + 0.0002])

    # Perform curve fitting
    try:
        popt, _ = curve_fit(gaussian_profile, wave, norm_flux, p0=initial_guess, bounds=bounds)
        b_opt, logN, lambda_0, z_opt = popt[1:5]  # Extract only the relevant parameters

        # Adjust components based on increments
        modified_components = [b_opt - increment_b, logN + increment_N, lambda_0, z_opt]
    except Exception as e:
        print(f"An error occurred during curve fitting: {e}")

    return modified_components, z_sys
