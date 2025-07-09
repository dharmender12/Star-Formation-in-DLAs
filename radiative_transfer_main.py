from radiative_transfer_mcmc import radiative_transfer

data_median = 'stacked_median_with_different_data.csv'
data_mean = 'weighted_mean_stacked_different_data.csv'
z_avg = 2.71
lo = 1205.
l1 = 1225.
radiative_transfer(data_median,data_mean,z_avg,lo,l1)