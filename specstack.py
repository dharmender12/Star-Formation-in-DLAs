##Python standard library
import sys
import os
import time
##python third party
import numpy
import scipy.stats as stats
from sklearn.utils import resample
import numpy as np
from scipy.stats import bootstrap

def restframe_normalised(spec, z, l0, l1, SNR, verbose):
    '''
    This function deredshift the wavelength grid
    and normalise the spectrum between l0 and l1

    Parameters
    ----------
    spec    str
            name of the spectrum to deredshift

    z       float
            redshift of the spectrum

    l0      float
            lower limit of the zone where we normalize

    l1      float
            upper limit ''' ''' ''' ''' ''' '''

    SNR     str
            SNR,l1_snr,l0_snr threshold of the SNR to compute in
            between snr_l0 and snr_l1

    verbose Bool
            True or false for printouts

    Returns
    -------
    wave_0  numpy array
            restframe wavelength
    flux_norm
            flux_density normalised in the normalisation region
    '''

    ##loadspectrum
    specdata = numpy.genfromtxt(spec,dtype='float').T
    wave = specdata[0]
    flux = specdata[1]

    ##deredshift the wavelength
    wave_0 = wave/(1+z) 

    ###SNR --> if the user gave something then we must compute it!
    if SNR != None:
        snr_threshold = float(SNR.split(',')[0])
        snr_l0 = float(SNR.split(',')[1])
        snr_lf = float(SNR.split(',')[2])
        reg_snr = numpy.where(numpy.logical_and(numpy.greater_equal(wave_0,snr_l0),\
                numpy.less_equal(wave_0,snr_lf)))
        
        ##compute snr
        mean = numpy.median(flux[reg_snr])
        std = numpy.std(flux[reg_snr])
        snr_data = mean/std

        if snr_data < snr_threshold:
            if verbose:
                print('The SNR is too low for spec %s (SNR=%s)'%(os.path.basename(spec), snr_data))
            return [], []


    if z>7:
        return [], []
    ###find the region where to normalize
    reg = numpy.where(numpy.logical_and(numpy.greater_equal(wave_0,l0),numpy.less_equal(wave_0,l1)))

    ###check if the region reports anything
    if len(reg[0]) == 0:
        if verbose:
            print('The region does not exist in spec %s (restframe), we skip it'%os.path.basename(spec))
        return [], []
    
    region = wave_0[reg]

    ###compute the median in the region
    med = numpy.mean(flux[reg])

    ##and devide all the spectrum by this amount
    # flux_norm = flux/med 
    flux_norm = flux
    return wave_0, flux_norm


def renorm(wave, flux,conti, z, l0, l1, verbose):
    '''
    This function deredshift the wavelength grid
    and normalise the spectrum between l0 and l1

    Parameters
    ----------
    wave    numpy array
            wavelength
    flux    numpy.array of flux

    z       float
            redshift of the spectrum
    l0      float
            lower limit of the zone where we normalize
    l1      float
            upper limit ''' ''' ''' ''' ''' '''

    Returns
    -------
    wave_0  numpy array
            restframe wavelength
    flux_norm
            flux_density normalised in the normalisation region
    '''

    ##deredshift the wavelength
    wave_0 = wave/(1+z) 
    flux = flux*(1+z)
    ###find the region where to normalize
    reg = numpy.where(numpy.logical_and(numpy.greater_equal(wave_0,l0),numpy.less_equal(wave_0,l1)))

    ###check if the region reports anything
    if len(reg[0]) == 0:
        if verbose:
            print('The region does not exist in spec %s(restframe), we skip'%spec)
        return [], []

    region = wave_0[reg]

    ###compute the median in the region
    med = numpy.median(flux[reg])

    ##and devide all the spectrum by this amount
    # flux_norm = flux/conti
    # flux_norm = flux/med 
    flux_norm = flux 
    # print(flux_norm)
    return wave_0, flux_norm



def regrid(restframe,grid,flux_err,tau,sigma):
    '''
    This function rebins all the spectra to the 
    new grid

    Parameters
    ----------
    restframe   list
                of all the spectra to stack

    grid        numpy.array
                new wave to interpolate the spectrum to

    Returns
    -------
    rebinned    list
                of all spectrum interpolated to the new grid
    '''
    print("Stacking is running please be patient :")
    
    clipped = []
    rebinned = []
    id_tau = numpy.where((tau >=10.)&(grid<1215.6701))[0]
    # id_tau = numpy.where((tau >=10.))[0]
    # id_tau = numpy.where((tau >=10)&(grid>=1215.6701))[0]
    for i in restframe:
        
        r = numpy.interp(grid, i[0], i[1])
        r_2 = numpy.interp(grid, i[0], i[1])
        
        r_1, low, upp = stats.sigmaclip(r_2, sigma, sigma)
        # print(r_1[id_tau][0])
        # n_samples =len(flux_err)
        # err =  np.sqrt(np.sum(flux_err**2)) /n_samples
        
        med = numpy.mean(r_1[id_tau])
        # weighted_mean = np.sum(med / (err ** 2)) / np.sum(1 / (err ** 2))
        weighted_mean = np.sum(med / (np.std(r_1) ** 2)) / np.sum(1 / (np.std(r_1) ** 2))
        # mean_r = np.nanmean(r_1)
        # # # if med > 1.0 : 
        r = r_2 - weighted_mean
        # print('the median is',med)
        # r = r_2 + med
        index_nan_min = numpy.where(grid<i[0][0])
        r[index_nan_min] = numpy.nan
        index_nan_max = numpy.where(grid>i[0][-1])
        r[index_nan_max] = numpy.nan
        # clipped.append(r_2[id_tau])
        rebinned.append(r)
        
        
    return rebinned, clipped




def stack(rebinned, sigma, flux_err):
    '''
    Function that stacks the spectra using a sigma-clipping method.

    Parameters
    -----------
    rebinned    list
                of rebinned spectra
                
    sigma       float
                sigma-clipping parameter
                
    flux_err    list
                list of flux errors for each spectrum

    Returns
    ------
    stacked     list
                stacked spectrum
                
    std         list
                associated standard deviation
    '''
    
    
    
    
    print("Bootstrapping is working")
    
    ## Convert to numpy array
    rebin = np.array(rebinned)
    specs = rebin.T
    
    err = []
    ### stacked_mean ######
    stacked_mean = []
    std_mean = []
    ermean = []
    N_mean = []
    #### stacked median ######
    stacked_median = []
    std_median = []
    ermedian = []
    N_median = []
    print('size of specs', len(specs))
    
    for i in specs:
        no_nan = np.count_nonzero(~np.isnan(i))
        index = np.where(np.isnan(i) == False)
        N_mean.append(no_nan)
        N_median.append(no_nan)
        
        if no_nan == 1:
            stacked_median.append(i[index][0])
            stacked_mean.append(i[index][0])
            std_median.append(i[index][0])
            std_mean.append(i[index][0])
            ermean.append(i[index][0])
            ermedian.append(i[index][0])
            err.append(i[index][0])
            
        elif no_nan < 10:
            stacked_median.append(np.nanmedian(i[index]))
            std_median.append(np.nanstd(i[index]))
            ermedian.append(np.nanstd(np.sqrt(no_nan)))
            stacked_mean.append(np.nanmean(i[index]))
            std_mean.append(np.nanstd(i[index]))
            ermean.append(np.nanstd(np.sqrt(no_nan)))
            err.append(np.nanstd(np.sqrt(no_nan)))
            
        else:
            c, low, upp = stats.sigmaclip(i[index], sigma, sigma)
            # print('the len of c', len(c), c.shape)

            num_bootstraps = 5000
            bootstrapped_means = []
            bootstrapped_median = []
            
            # res = bootstrap((c,), np.std, confidence_level=0.9, n_resamples=num_bootstraps, random_state=np.random.default_rng())
            # ci_l, ci_u = res.confidence_interval
            # Extracting bootstrap distribution from bootstrap samples
            # bootstrap_distribution = res.bootstrap_distribution

            # resampled_mean = np.nanmean(c, axis=0)
            # resampled_median = numpy.nanmedian(c, axis=0) 
            # bootstrapped_median.append(resampled_median)       
            #####Append the standard error from the bootstrap result
            # err.append(res.standard_error)
            for _ in range(num_bootstraps):
                resampled_flux = resample(c, replace=True)
                # resampled_error = resample(c, replace=True)
                # Compute the median for each resampled dataset
                resampled_mean = np.nanmean(resampled_flux, axis=0)
                resampled_median = numpy.nanmedian(resampled_flux, axis=0)
                # Append the median to the list of bootstrapped medians
                bootstrapped_median.append(resampled_median)
                
                
                
                # flux_boot = bootstrap(c, statistic=np.median, 
                #          n_resamples=1000, confidence_level=0.68, method='percentile')
                # print(bootstrapped_medians)
                ##### --------------- weighted mean --------------- ####################### 
                std_err = np.nanstd(resampled_flux,axis=0)
                # weighted_mean = np.sum(resampled_mean / (err ** 2)) / np.sum(1 / (std_err ** 2))
                weighted_mean = np.sum(resampled_mean / (std_err ** 2)) / np.sum(1 / (std_err ** 2))
                    # err.append(std_err)
                    # Calculate weighted mean
                    # print(weighted_mean)
                bootstrapped_means.append(weighted_mean)
            # res = bootstrap(c, np.std, confidence_level=0.9,axis=0,n_resamples=1000,
            #     random_state= np.random.default_rng())
            ########## median ##############
            # stacked.append(np.nanmedian(c,axis=0))
            # stdeviation.append(np.nanstd(c,axis=0))
            #### Bootstrapped Median ##############
            stacked_median.append(numpy.nanmedian(bootstrapped_median))
            std_median.append(numpy.nanstd(bootstrapped_median,axis = 0))
            ermedian.append(numpy.nanstd(bootstrapped_median)/numpy.count_nonzero(~numpy.isnan(bootstrapped_median)))
                # print(numpy.nanmedian(bootstrapped_medians).shape)
            ## Bootstrapped Mean ##############
            stacked_mean.append(np.nanmean(bootstrapped_means))
            std_mean.append(np.nanstd(bootstrapped_means))
            # std_mean.append(std_err)
            ermean.append(np.nanstd(bootstrapped_means) / np.count_nonzero(~np.isnan(bootstrapped_means)))
            ### ------- Root mean square error ----------###############
            # Calculate flux error
            
            
            # print('The sample size is :',n_samples,flux_err)
            # res_flux  = (resampled_flux,axis = 0)
            # flux_err = np.sqrt(np.sum((resampled_flux - resampled_mean) ** 2) / n_samples) * 1 / np.sqrt(n_samples)
            
            # print(er)
            n_samples =len(flux_err)
            flux_er =  np.sqrt(np.sum(flux_err**2)) /n_samples
            err.append(flux_er)
        

    return stacked_mean, std_mean, ermean, N_mean, stacked_median, std_median, ermedian, N_median, err


# def stack(rebinned, sigma,flux_err):
#     '''
#     Function that stacks the spectra
#     using a sigma-clipping method.

#     Parameters
#     -----------
#     rebinned    list
#                 of rebinned spectra
                
#     sigma       float
#                 sigma-clipping parameter
#     return
#     ------
#     stacked     list
#                 stacked spectrum
#     std         list
#                 associated standard deviation
#     '''
#     import numpy as np
#     from scipy import stats
#     from sklearn.utils import resample
#     from scipy.stats import bootstrap
    
#     print("Bootstrapping is working ")
#     ## Convert to numpy array
#     rebin = np.array(rebinned)
#     specs = rebin.T
#     # print('the spectrum values are :',specs)
#     err = []
# 	### stacked_mean ######
#     stacked_mean = []
#     std_mean = []
#     ermean = []
#     N_mean = []
#     #### stacked median ######
#     stacked_median = []
#     std_median = []
#     ermedian = []
#     N_median = []
#     print('size of specs',len(specs))
#     for i in specs:
#         no_nan = np.count_nonzero(~np.isnan(i))
#         index = np.where(np.isnan(i) == False)
#         N_mean.append(no_nan)
#         N_median.append(no_nan)
#         if no_nan == 1:
#             stacked_median.append(i[index][0])
#             stacked_mean.append(i[index][0])
#             std_median.append(i[index][0])
#             std_mean.append(i[index][0])
#             ermean.append(i[index][0])
#             ermedian.append(i[index][0])
#             err.append(i[index][0])
#         elif no_nan < 10:
#             stacked_median.append(np.nanmedian(i[index]))
#             std_median.append(np.nanstd(i[index]))
#             ermedian.append(np.nanstd(np.sqrt(no_nan)))
#             stacked_mean.append(np.nanmean(i[index]))
#             std_mean.append(np.nanstd(i[index]))
#             ermean.append(np.nanstd(np.sqrt(no_nan)))
#             err.append(np.nanstd(np.sqrt(no_nan)))
#         else:
#             c, low, upp = stats.sigmaclip(i[index], sigma, sigma)
#             print('the len of c',len(c),c.shape)

#             num_bootstraps = 5000
#             bootstrapped_means = []
#             bootstrapped_median = []
        
 
#             res = bootstrap((c,), np.median, confidence_level=0.9, n_resamples=num_bootstraps, random_state=np.random.default_rng())
#             bootstrap_distribution = res.bootstrap_distribution

#             # Compute the mean and median of the bootstrap distribution
#             bootstrap_mean = np.mean(bootstrap_distribution)
#             bootstrap_median = np.median(bootstrap_distribution)
            
#             bootstrapped_means.append(bootstrap_mean)
#             bootstrapped_median.append(bootstrap_median)

#             # Append the standard error from the bootstrap result
#             err.append(res.standard_error)

#             print(res, res.confidence_interval)

   
            # for _ in range(num_bootstraps):
            #     resampled_flux = resample(c, replace=True)
            #     # resampled_error = resample(c, replace=True)
            #     # Compute the median for each resampled dataset
            #     resampled_mean = np.nanmean(resampled_flux, axis=0)
            #     resampled_median = numpy.nanmedian(resampled_flux, axis=0)
            #     # Append the median to the list of bootstrapped medians
            #     bootstrapped_median.append(resampled_median)
                
                
                
            #     # flux_boot = bootstrap(c, statistic=np.median, 
            #     #          n_resamples=1000, confidence_level=0.68, method='percentile')
            #     # print(bootstrapped_medians)
            # ##### --------------- weighted mean --------------- ####################### 
            #     std_err = np.nanstd(resampled_flux,axis=0)
            #     weighted_mean = np.sum(resampled_mean / (std_err ** 2)) / np.sum(1 / (std_err ** 2))
            #     # err.append(std_err)
            #     # Calculate weighted mean
            #     # print(weighted_mean)
            #     bootstrapped_means.append(weighted_mean)
            # # res = bootstrap(c, np.std, confidence_level=0.9,axis=0,n_resamples=1000,
            # #     random_state= np.random.default_rng())
            # ########## median ##############
            # # stacked.append(np.nanmedian(c,axis=0))
            # # stdeviation.append(np.nanstd(c,axis=0))
            # #### Bootstrapped Median ##############
            # stacked_median.append(numpy.nanmedian(bootstrapped_median))
            # std_median.append(numpy.nanstd(bootstrapped_median,axis = 0))
            # ermedian.append(numpy.nanstd(bootstrapped_median)/numpy.count_nonzero(~numpy.isnan(bootstrapped_median)))
            #     # print(numpy.nanmedian(bootstrapped_medians).shape)
            # ## Bootstrapped Mean ##############
            # stacked_mean.append(np.nanmean(bootstrapped_means))
            # std_mean.append(np.nanstd(bootstrapped_means))
            # # std_mean.append(std_err)
            # ermean.append(np.nanstd(bootstrapped_means) / np.count_nonzero(~np.isnan(bootstrapped_means)))
            # ### ------- Root mean square error ----------###############
            # # Calculate flux error
            # n_samples =len(flux_err)
            
            # # print('The sample size is :',n_samples,flux_err)
            # # res_flux  = (resampled_flux,axis = 0)
            # # flux_err = np.sqrt(np.sum((resampled_flux - resampled_mean) ** 2) / n_samples) * 1 / np.sqrt(n_samples)
            # er =  np.sqrt(np.sum(flux_err**2)) /n_samples
            # # print(er)
            
            # err.append(res.standard_error)
      

    # return np.array(stacked_mean), np.array(std_mean), np.array(ermean), np.array(N_mean), np.array(stacked_median),np.array(std_median),np.array(ermedian),np.array(N_median),np.array(err)



