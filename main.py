import os
import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from line_choice import velocity_plot  # Ensure this is the correct import path for your velocity_plot function
from process_file import process_file
from scipy.optimize import curve_fit

from matplotlib import rcParams
# Set the default font family
rcParams['font.family'] = 'DejaVu Sans'
# Define the filtering conditions
value_1 = 20.8
value_2 = 4.0
res = 2000
# Prepare to store results
results = []

file = '/home/dharmender/New_dla/dla_with_non_matching_ids.csv'
# file = '/home/dharmender/dla_new_file/DLA/coadded_specfname_dla_with_cnr.csv'
# Load the DLA data
dla_data_with_cnr = pd.read_csv(file)
filepath = '/home/dharmender/New_dla/Coadd_Spec/'
# filepath = '/home/dharmender/dla_new_file/DLA/Coadd_Spec/'
# Define the filtering conditions
filtered_data = dla_data_with_cnr
# filtered_data = dla_data_with_cnr[(dla_data_with_cnr['NHIfit'] >= value_1) & (dla_data_with_cnr['CNR'] >= value_2)]
print('the length of data is : ',len(filtered_data))
# Initialize counter
counter = 0

# Loop through the filtered DLA data  

for i_file in range(len(filtered_data)):
# for i_file in range(528,len(filtered_data)):
# for i_file in range(760,len(filtered_data)):
    ra = filtered_data['Radeg'].iloc[i_file]
    dec = filtered_data['Dedeg'].iloc[i_file]
    if filtered_data['coadd_bossname'].iloc[i_file] is not None and not pd.isnull(filtered_data['coadd_bossname'].iloc[i_file]) and filtered_data['coadd_bossname'].iloc[i_file] != '0':
        filename = f"/home/dharmender/New_dla/Coadd_Spec/{filtered_data['coadd_bossname'].iloc[i_file]}"
        z_dla = filtered_data['zCNN'].iloc[i_file]
    elif filtered_data['coadd_sdssname'].iloc[i_file] is not None and not pd.isnull(filtered_data['coadd_sdssname'].iloc[i_file]) and filtered_data['coadd_sdssname'].iloc[i_file] != '0':
        filename = f"/home/dharmender/New_dla/Coadd_Spec/{filtered_data['coadd_sdssname'].iloc[i_file]}"
        z_dla = filtered_data['zCNN'].iloc[i_file]
    else:
        filename = f"/home/dharmender/New_dla/DLASTACk/{filtered_data['Generated_Specfname'].iloc[i_file]}"
        z_dla = filtered_data['zCNN'].iloc[i_file]

    z_cnn = filtered_data['zCNN'].iloc[i_file]
    id = filtered_data['ID_1'].iloc[i_file]
    logNHI = filtered_data['NHIfit'].iloc[i_file], 0.1
    nh1 = filtered_data['NHIfit'].iloc[i_file]
    cnr = filtered_data['CNR'].iloc[i_file]
    # Read FITS file
    with fits.open(filename) as hdul:
        data = hdul[1].data  # Assuming the wavelength is in the first HDU
        wave = data['WAVE'][0]
        flux = data['FLUX'][0]  # Assuming STACK_MEDIAN contains optical depth
        err = data['ERROR'][0]
        conti_tm = data['CONTI_TM'][0]

    normalised_flux = flux / conti_tm
    normalised_flux[np.isinf(normalised_flux)] = 0  # Replace inf with 0
    
    # Get user choice for spectrum plot

    norm_err = (err / conti_tm)
    norm_err[np.isinf(norm_err)] = 0  # Replace inf with 0
    
    # Call the process_file function with your parameters
    dataset,z_correct, z_sys, primary_lines, centroid_value  = process_file(filename, id, wave, normalised_flux, z_cnn, res, norm_err, logNHI, increment_b=1.0, increment_N=0.5, modify_params=True, max_attempts=10)
    
    
    # Save the result to the list
    results.append([id, ra,dec,os.path.basename(filename), z_cnn, z_correct, z_sys, primary_lines, centroid_value,cnr,nh1])
    
    # Save results to CSV after processing each file
    results_df = pd.DataFrame(results, columns=['ID','Radeg','Dedeg', 'Filename', 'zCNN', 'z_choice', 'z_metal', 'metal', 'Line_center','CNR','logNHI'])
    results_df.to_csv('filtered_dla_data_cnr_4.csv', index=False)
    
    print("File number and :",counter+1,id,os.path.basename(filename),z_correct)
    # Increment counter
    counter += 1
