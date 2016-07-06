#!~/bin/python
# Python spectral indicies
# Corey Mutnik 5/31/16

# in ipython type history to see what you typed in that session
#history

from pysynphot import spectrum, observation
import numpy.polynomial.polynomial as poly
import matplotlib.gridspec as gridspec
import scipy.optimize as optimized
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits
from glob import glob
import numpy as np
import pysynphot
import os


def rebin_spec(wave, specin, wavenew):
    '''Rebin data for direct comparison'''
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux = specin)
    f = np.ones(len(wave)) # returns array of ones with size given
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavenew, force='taper')
    return obs  
def normalize_flux(wave_array, flux_array):
    micron_value = np.where(abs(wave_array - 2.20) == min(abs(wave_array - 2.20)))[0][0]
    norm_den = (float)((wave_array[micron_value] * flux_array[micron_value])**(-1))
    norm_flux = []
    for i in range(0, len(wave_array)):
        norm_flux.append(wave_array[i] * flux_array[i] * norm_den)
    return np.array(norm_flux)  
def normalize_spt_feature(wave_array, flux_array, normval):
    micron_value = np.where(abs(wave_array - normval) == min(abs(wave_array - normval)))[0][0]
    norm_den = (float)((wave_array[micron_value] * flux_array[micron_value])**(-1))
    norm_flux = []
    for i in range(0, len(wave_array)):
        norm_flux.append(wave_array[i] * flux_array[i] * norm_den)
    return np.array(norm_flux) 
def eval_poly(x, a,b,c,d,e):
    return (a + b*x + c*x**2 + d*x**3 + e*x**4)
def eval_linfit(x, m,b):
    return(m*x + b)

# degree of polynomial
polydeg = 4 

'''
def Ca_II_EW_has_clunky(inwave,influx):
    global EW_caii
    # Ca ii (0.866 \mum)
    feature_lamb=0.866
    # Feature limits: 0.860-0.875
    window = (0.860 <= inwave) & (0.875 >= inwave)

    # mask wavelength and flux, so only values in window remain
    wave_caii = inwave[window]
    flux_caii = influx[window]

    ###
    # NORMALIZE TO FEATURE???
    ##
    # ANS: NO
    ###--------------------------------------------------------------------------------------------------------------
    #flux_norm2feat = normalize_spt_feature(wave_array=wave_caii, flux_array=flux_caii, normval=feature_lamb)
    #----------------------------------------------------------------------------------------------------------------

    # First continuum level: 0.862-0.864
    continuum1 = (0.862 <= inwave) & (0.864 >= inwave)
    #Second continuum level: 0.870-0.873
    continuum2 = (0.870 <= inwave) & (0.873 >= inwave)

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    


    # only use region between first and second continuum regions
    feature_within_limits = (0.864 <= inwave) & (0.870 >= inwave)
    wave_caii_within_limits = inwave[feature_within_limits]
    flux_caii_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_caii_cont = eval_linfit(wave_caii_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_caii_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_caii = (  (1.0 - (flux_caii_within_limits / flux_caii_cont)) * delta_wave).sum()

    def plot_CaII_window():
        plt.clf()
        plt.plot(wave_caii,flux_caii, 'ob')
        plt.plot(wave_caii,flux_caii, 'b')
        plt.plot(wave4linfit,flux4linfit,'ok')
        plt.plot(inwave[feature_within_limits],influx[feature_within_limits], 'og')
        plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
        plt.savefig(standard_file+'_R09_CaII_window_fit.png')
        #plt.show()
    plot_CaII_window()
    def EW_clunky_method():
        EW_clunky = 0.
        for jj in range(len(flux_caii)):
            f_lambi = flux_caii[jj]
            fc_lambi = eval_linfit(wave_caii[jj], linparam[0], linparam[1])
            if jj==0:
                # no point before the first one, so double the distance
                # between the first and second points
                #print 'first:', jj
                dellambi = (wave_caii[1] - wave_caii[0]) * 2.
            elif jj==(len(flux_caii)-1):
                # no point after the last one, so double the distance
                # between the last and second to last point
                #print 'last:', jj
                dellambi = (wave_caii[-1] - wave_caii[-2]) * 2.
            else:
                # \Delta\lambda_i = \lambda{i+1} - \lambda{i-1}
                dellambi = wave_caii[jj+1] - wave_caii[jj-1]
            dellambi*=0.5*1e4
            #dellambi*=0.5
            EW_clunky += (1. - (f_lambi/fc_lambi))*dellambi
        #print 'myold_clunky_EW_clunky_method() EW_clunky value:', EW_clunky
        #>> EW_clunky value: 0.0121197321471
        global EW_clunky
    #EW_clunky_method()
    def integrate_gaussian_unfinished():
        # Gaussian fit
        def gaussian(x, a, b, c, d):
            #a = (height of the peak - continuum), b = peak center, c = width of peak, d = continuum
            y = a * np.exp((-(x-b)**2)/(2*c**2)) + d
            return y

        #a=5e-11 - 6e-11
        a=5*10**-11 - 6*10**-11
        b,c,d=0.8665, .005, 6.5*10**-11

        param, covar = optimized.curve_fit(gaussian, wave_caii, flux_caii, p0=[a,b,c,d])
        #print param 

        y_plot = gaussian(wave_caii, param[0], param[1], param[2], param[3]) # plot fit against the x_data
        gaussians, = plt.plot(wave_caii, y_plot,'.r', label='Gaussian Fits')

        from scipy.integrate import quad
        I = quad(gaussian, min(wave_caii), max(wave_caii), args=(param[0], param[1], param[2], param[3]))
    def from_jessica():
        global EW_caii
        flux_caii_cont = eval_linfit(wave_caii, linparam[0], linparam[1])
        delta_wave = np.diff(wave_caii)
        delta_wave = np.append(delta_wave, delta_wave[-1])
        EW = (  (1.0 - (flux_caii / flux_caii_cont)) * delta_wave).sum()
        #print 'EW value:', EW
        #print EW
        #print flux_caii_cont[30]
        #print flux_caii[30]
        #(flux_caii/flux_caii_cont)[30]
        
        # convert microns to Angs
        delta_wave *= 1e4
        EW = (  (1.0 - (flux_caii / flux_caii_cont)) * delta_wave).sum()
        #print EW
        #>> 4.8754441186

        def quickplot():
            plt.clf()
            plt.plot(wave_caii, flux_caii/flux_caii_cont, 'k.')
            plt.title('$f(\lambda) / f_c(\lambda)$')
            plt.savefig(standard_file+'_CaII_jess_quickplot.png')
        quickplot() 

        def quickplot_flipped():
            plt.clf()
            plt.plot(wave_caii, 1.0 - (flux_caii/flux_caii_cont), 'k.')
            plt.title('1. - ($f(\lambda) / f_c(\lambda)$)')
            plt.savefig(standard_file+'_CaII_jess_quickplotflipped.png')
        quickplot_flipped() 
    
    

        # check header
        #head = fits.getheader(standard_dir)
        #print head
        
        def only_sum_between_limits():
            EW_tmp = (  (1.0 - (flux_caii / flux_caii_cont)) * delta_wave)
            #EW_tmp[15:45].sum()
            #idx = np.where((wave_caii >= 0.864) & (wave_caii <= 0.870))[0]
            #EW_tmp[idx].sum()
            EW_tmp[15:80].sum()
            idx = np.where((wave_caii > 0.864) & (wave_caii < 0.870))[0]
            EW_tmp[idx].sum()
            #idx = np.where((wave_caii => 0.864) & (wave_caii <= 0.870))[0]
            #EW_tmp[idx].sum()
            #idx
            #len(EW)
            #len(EW_tmp)

        flux_caii_cont = eval_linfit(wave_caii, linparam[0], linparam[1])
        delta_wave = np.diff(wave_caii)
        delta_wave = np.append(delta_wave, delta_wave[-1])
        EW = (  (1.0 - (flux_caii / flux_caii_cont)) * delta_wave).sum()
        #print EW
        # convert wavelength from microns to angstrom
        delta_wave *= 1e4
        EW = (  (1.0 - (flux_caii / flux_caii_cont)) * delta_wave).sum()
        #print EW


        
        ##EW_tmp = (  (1.0 - (flux_caii / flux_caii_cont)) * delta_wave)
        ###EW_tmp[15:45].sum()
        ###idx = np.where((wave_caii >= 0.864) & (wave_caii <= 0.870))[0]
        ###EW_tmp[idx].sum()
        ##EW_tmp[15:80].sum()
        ##idx = np.where((wave_caii > 0.864) & (wave_caii < 0.870))[0]
        ##EW_tmp[idx].sum()
        ###idx = np.where((wave_caii >= 0.864) & (wave_caii <= 0.870))[0]
        ###EW_tmp[idx].sum()
        ##idx
        ###len(EW)
        ###len(EW_tmp)
        ##
        ##def test_see_where_works_same_as_my_method():
        ##    idx = np.where((wave_caii > 0.864) & (wave_caii < 0.870))[0]
        ##    plt.clf()
        ##    plt.plot(wave_caii,flux_caii, '.b')
        ##    plt.plot(wave4linfit,flux4linfit,'.k')
        ##    plt.plot(inwave[feature_within_limits],influx[feature_within_limits], 'og')
        ##    plt.plot(wave_caii[idx],flux_caii[idx], '.r')
        ##    #plt.show()
        ###test_see_where_works_same_as_my_method()
        ###global EW, EW_tmp
        
        EW_caii = EW
    #from_jessica()



    #print '\n\nCa II\nclunky EW value:\t', EW_clunky , '\nrefined EW (Angstrom):\t', EW#, '\nEW_tmp:\t', EW_tmp[idx].sum()
# 160609 - has issues with variance...output nan
def Ca_II_EW(inwave,influx,seqnum='_'):
    global EW_caii, variance_caii
    # Ca II (0.866 \mum)
    # Feature limits: 0.860-0.875
    low_lim,high_lim = 0.860, 0.875
    # First continuum level: 0.862-0.864
    fc1,fc2 = 0.862, 0.864
    # Second continuum level: 0.870-0.873
    sc1,sc2 = 0.870, 0.873

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # change here -----------------------------------------------------------------------------------------------------------
    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_caii = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5
    variance_caii = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    var = np.sqrt(variance_caii)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_caii, var
    #------------------------------------------------------------------------------------------------------------------------
def Mg_I_1tac7_EW(inwave,influx,seqnum='_'):
    global EW_mgi_1tac7, variance_mgi_1tac7
    # Mg i (1.711 \mum)
    # Feature limits: 1.695-1.726
    low_lim,high_lim = 1.695, 1.726
    # First continuum level: 1.702-1.708
    fc1,fc2 = 1.702, 1.708
    #Second continuum level: 1.715-1.720
    sc1,sc2 = 1.715, 1.720

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_mgi_1tac7 = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    #--------------- ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5
    variance_mgi_1tac7 = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    var = np.sqrt(variance_mgi_1tac7)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_mgi_1tac7, var
def Mg_I_1tac5_EW(inwave,influx,seqnum='_'):
    global EW_mgi_1tac5, variance_mgi_1tac5
    # Mg i (1.485 \mum)
    # Feature limits: 1.475-1.4975
    low_lim,high_lim = 1.475, 1.4975
    # First continuum level: 1.4775-1.485
    fc1,fc2 = 1.4775, 1.485
    #Second continuum level: 1.491-1.497
    sc1,sc2 = 1.491, 1.497

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # change here -----------------------------------------------------------------------------------------------------------
    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_mgi_1tac5 = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5
    variance_mgi_1tac5 = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    var = np.sqrt(variance_mgi_1tac5)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_mgi_1tac5, var
    #------------------------------------------------------------------------------------------------------------------------
def Al_I_EW(inwave,influx,seqnum='_'):
    global EW_ali, variance_ali
    # Al i (1.313 \mum)
    # Feature limits: 1.300-1.330
    low_lim,high_lim = 1.300,1.330
    # First continuum level: 1.305-1.309
    fc1,fc2 = 1.305,1.309
    # Second continuum level: 1.320-1.325
    sc1,sc2 = 1.320,1.325

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # change here -----------------------------------------------------------------------------------------------------------
    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_ali = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5
    variance_ali = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    var = np.sqrt(variance_ali)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_ali, var
    #------------------------------------------------------------------------------------------------------------------------
def Na_I_1tac14_EW(inwave,influx,seqnum='_'):
    global EW_nai, variance_nai
    # Na i (1.14 \mum)
    # Feature limits: 1.120-1.160
    low_lim,high_lim = 1.120,1.160
    # First continuum level: 1.125-1.130
    fc1,fc2 = 1.125,1.130
    #Second continuum level: 1.150-1.160
    sc1,sc2 = 1.150,1.160

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # change here -----------------------------------------------------------------------------------------------------------
    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_nai = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5
    variance_nai = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    var = np.sqrt(variance_nai)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_nai, var
    #------------------------------------------------------------------------------------------------------------------------
def Na_I_2tac2_EW(inwave,influx):
    global EW_nai_2tac2, variance_nai_2tac2
    # Na i (2.206 \mum) 
    # Feature limits: 2.185-2.230
    low_lim, high_lim = 2.185,2.230
    # First continuum level: 2.192-2.198
    fc1,fc2 = 2.192,2.198
    # Second continuum level: 2.213-2.220
    sc1,sc2 = 2.213,2.220

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # change here -----------------------------------------------------------------------------------------------------------
    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_nai_2tac2 = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5
    variance_nai_2tac2 = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    var = np.sqrt(variance_nai_2tac2)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_nai_2tac2, var
    #------------------------------------------------------------------------------------------------------------------------
'''

# 160616 - fixed variance issues
def Ca_II_EW(inwave,influx):
    global EW_caii, variance_caii
    # Ca II (0.866 \mum)
    # Feature limits: 0.860-0.875
    low_lim,high_lim = 0.860, 0.875
    # First continuum level: 0.862-0.864
    fc1,fc2 = 0.862, 0.864
    # Second continuum level: 0.870-0.873
    sc1,sc2 = 0.870, 0.873

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_caii = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5

    #variance_caii = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    # variance above fails (gives nan) due to flux_cont^-4, so modify below
    variance_caii = (( (sig / flux_cont)**2 + ((sigc * flux_within_limits / flux_cont)**2 * flux_cont**-2) ) * delta_wave**2).sum()
    var = np.sqrt(variance_caii)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_caii, var
def Mg_I_1tac7_EW(inwave,influx):
    global EW_mgi_1tac7, variance_mgi_1tac7
    # Mg i (1.711 \mum)
    # Feature limits: 1.695-1.726
    low_lim,high_lim = 1.695, 1.726
    # First continuum level: 1.702-1.708
    fc1,fc2 = 1.702, 1.708
    #Second continuum level: 1.715-1.720
    sc1,sc2 = 1.715, 1.720

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_mgi_1tac7 = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    #--------------- ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5

    #variance_mgi_1tac7 = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    # variance above fails (gives nan) due to flux_cont^-4, so modify below
    variance_mgi_1tac7 = (( (sig / flux_cont)**2 + ((sigc * flux_within_limits / flux_cont)**2 * flux_cont**-2) ) * delta_wave**2).sum()

    var = np.sqrt(variance_mgi_1tac7)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_mgi_1tac7, var
def Mg_I_1tac5_EW(inwave,influx):
    global EW_mgi_1tac5, variance_mgi_1tac5
    # Mg i (1.485 \mum)
    # Feature limits: 1.475-1.4975
    low_lim,high_lim = 1.475, 1.4975
    # First continuum level: 1.4775-1.485
    fc1,fc2 = 1.4775, 1.485
    #Second continuum level: 1.491-1.497
    sc1,sc2 = 1.491, 1.497

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_mgi_1tac5 = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5

    #variance_mgi_1tac5 = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    # variance above fails (gives nan) due to flux_cont^-4, so modify below
    variance_mgi_1tac5 = (( (sig / flux_cont)**2 + ((sigc * flux_within_limits / flux_cont)**2 * flux_cont**-2) ) * delta_wave**2).sum()

    var = np.sqrt(variance_mgi_1tac5)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_mgi_1tac5, var
def Al_I_EW(inwave,influx):
    global EW_ali, variance_ali
    # Al i (1.313 \mum)
    # Feature limits: 1.300-1.330
    low_lim,high_lim = 1.300,1.330
    # First continuum level: 1.305-1.309
    fc1,fc2 = 1.305,1.309
    # Second continuum level: 1.320-1.325
    sc1,sc2 = 1.320,1.325

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_ali = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5

    #variance_ali = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    # variance above fails (gives nan) due to flux_cont^-4, so modify below
    variance_ali =(( (sig / flux_cont)**2 + ((sigc * flux_within_limits / flux_cont)**2 * flux_cont**-2) ) * delta_wave**2).sum()

    var = np.sqrt(variance_ali)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_ali, var
def Na_I_1tac14_EW(inwave,influx):
    global EW_nai, variance_nai
    # Na i (1.14 \mum)
    # Published Feature limits: 1.120-1.160
    low_lim,high_lim = 1.120,1.160
    # First continuum level: 1.125-1.130
    fc1,fc2 = 1.125,1.130
    #Second continuum level: 1.150-1.160
    sc1,sc2 = 1.150,1.160

    ##
    # Corrected Feature Limits
    ##
    FL1, FL2 = 1.137,1.1428

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (FL1 <= inwave) & (FL2 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_nai = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5

    #variance_nai = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    # variance above fails (gives nan) due to flux_cont^-4, so modify below
    variance_nai =(( (sig / flux_cont)**2 + ((sigc * flux_within_limits / flux_cont)**2 * flux_cont**-2) ) * delta_wave**2).sum()

    var = np.sqrt(variance_nai)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_nai, var
def Na_I_2tac2_EW(inwave,influx):
    global EW_nai_2tac2, variance_nai_2tac2
    # Na i (2.206 \mum) 
    # Feature limits: 2.185-2.230
    low_lim, high_lim = 2.185,2.230
    # First continuum level: 2.192-2.198
    fc1,fc2 = 2.192,2.198
    # Second continuum level: 2.213-2.220
    sc1,sc2 = 2.213,2.220

    # define window of region (was used for plotting)
    #window = (low_lim <= inwave) & (high_lim >= inwave)
    #wave_window = inwave[window]
    #flux_window = influx[window]

    # continuum regions
    continuum1 = (fc1 <= inwave) & (fc2 >= inwave)
    continuum2 = np.where((sc1 <= inwave) & (sc2 >= inwave))[0]

    # mask regions for unweighted linear fit
    wave4linfit = np.append(inwave[continuum1], inwave[continuum2])
    flux4linfit = np.append(influx[continuum1], influx[continuum2])

    # initial guesses
    a,b = 1.,1.
    
    linparam,lincovar = optimized.curve_fit(eval_linfit, wave4linfit, flux4linfit, p0=[a,b])
    #print linparam
    
    lifit_plot = eval_linfit(wave4linfit, linparam[0], linparam[1]) # plot fit against the x_data
    #limits_continuum, = plt.plot(wave4linfit, lifit_plot, label='Linear Fit', color='red')
    
    # only use region between first and second continuum regions
    feature_within_limits = (fc2 <= inwave) & (sc1 >= inwave)
    wave_within_limits = inwave[feature_within_limits]
    flux_within_limits = influx[feature_within_limits]

    # calculate EW
    flux_cont = eval_linfit(wave_within_limits, linparam[0], linparam[1])
    delta_wave = np.diff(wave_within_limits)
    # add the last value to the end of array, so they are same shape
    delta_wave = np.append(delta_wave, delta_wave[-1])

    # convert wavelength from microns to angstrom
    delta_wave *= 1e4
    EW_nai_2tac2 = (  (1.0 - (flux_within_limits / flux_cont)) * delta_wave).sum()

    # ERROR CALCULATION
    sig = np.std(flux4linfit)
    sigc = sig * (lincovar[0,0]**2 + lincovar[1,1]**2)**0.5

    #variance_nai_2tac2 = (( (sig**2 / flux_cont**2) + (flux_within_limits**2 / flux_cont**4)*sigc**2 ) * delta_wave**2).sum()
    # variance above fails (gives nan) due to flux_cont^-4, so modify below
    variance_nai_2tac2 = (( (sig / flux_cont)**2 + ((sigc * flux_within_limits / flux_cont)**2 * flux_cont**-2) ) * delta_wave**2).sum()
    
    var = np.sqrt(variance_nai_2tac2)
    print 'sig   sigc   variance   var'
    print sig, sigc, variance_nai_2tac2, var

# make list of all fits files
#dir_path = u'/Users/cmutnik/work/astro/plots/compare/z_rayner_all'
dir_path = '/Users/cmutnik/work/astro/plots/compare/z_rayner_all'
globpath = os.path.join(dir_path, '*.fits')
filelist = glob(globpath)
filelist.sort()

# set file to write EW values to
outfile = open('EW_of_rayner_stars_fixed_NaI_1tac14.txt', 'w')

# write header to file, for colnames
#outfile.write('#standard_file, EW_caii, EW_nai, EW_nai_2tac2, EW_ali, EW_mgi_1tac5, EW_mgi_1tac7, variance_mgi_1tac7\n')
outfile.write('#standard_file, EW_caii, variance_caii, EW_nai, variance_nai, EW_ali, variance_ali, EW_mgi_1tac5, variance_mgi_1tac5, EW_mgi_1tac7, variance_mgi_1tac7, EW_nai_2tac2, variance_nai_2tac2\n')
for loopy in range(len(filelist)):
    standard_file = os.path.basename(filelist[loopy])[:-5]
    # open IRTF data
    data_irtf = fits.getdata(filelist[loopy])
    wave_irtf = data_irtf[0]
    flux_irtf = data_irtf[1]

    flux_mask_irtf = flux_irtf[np.logical_not(np.isnan(flux_irtf))]
    wave_mask_irtf = wave_irtf[np.logical_not(np.isnan(flux_irtf))]

    # get EW values for each spectral line
    Ca_II_EW(wave_mask_irtf,flux_mask_irtf)
    Na_I_1tac14_EW(wave_mask_irtf,flux_mask_irtf)
    Al_I_EW(wave_mask_irtf,flux_mask_irtf)
    Na_I_2tac2_EW(wave_mask_irtf,flux_mask_irtf)
    Mg_I_1tac5_EW(wave_mask_irtf,flux_mask_irtf)
    Mg_I_1tac7_EW(wave_mask_irtf,flux_mask_irtf)

    # write values to file
    #string_um_up = (standard_file,EW_caii,EW_nai,EW_nai_2tac2,EW_ali,EW_mgi_1tac5,EW_mgi_1tac7,variance_mgi_1tac7)
    string_um_up = (standard_file,EW_caii,variance_caii,EW_nai,variance_nai,EW_ali,variance_ali,EW_mgi_1tac5,variance_mgi_1tac5,EW_mgi_1tac7,variance_mgi_1tac7,EW_nai_2tac2,variance_nai_2tac2)
    string_um_up = str(string_um_up)+'\n'
    outfile.write(string_um_up)
    # print where im at in loop
    print 'Finished: ', (loopy+1),'/',(len(filelist))
outfile.close()
