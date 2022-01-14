import xarray as xr
import numpy as np
from scipy.stats.distributions import chi2
from scipy import signal
import xrft
import matplotlib.pyplot as plt

def spec_xr(da, dim, window=True):
    """ Computes power spectrum of an DatArray using xrft
    
    Parameters
    ----------
    da : `xarray.DataArray`
    dim : str
        The dimensions along which to take the transformation
    window : boolean
        See xrft documentation
        
    Returns
    -------
    da : `xarray.DataArray`
        power spectrum of the input DataArray
    
    """
    atrack = da[dim] - da[dim].min()
    da[dim] = atrack
    tmp = da.interpolate_na(dim=dim)
    sp = xrft.power_spectrum(
        tmp.dropna(dim=dim),dim=dim, detrend='linear',
        window=window).compute()
    freq_dim = 'freq_{}'.format(dim)
    ind = sp[freq_dim] > 0
    return 2*sp.isel(**{freq_dim: ind})


def daniell(da, dim, npoints, start_freq=None):
    """ Averages the power spectrum using the Daniell method
   
   Parameters
    ----------
    da : `xarray.DataArray`
        power spectrum from xrft
    dim : str
        The dimensions along which to take the transformation
    npoints : int
        number of points to average
    start_freq: float, optional
        performs the averaging only for points above start_freq. If None,
        it will perform the averaging for all frequencies
    
    Returns
    -------
    spec_dani : `xarray.DataArray`
        smoothed power spectrum
    box_width : `xarray.DataArray`
        number of points in the original spectrum that were 
        avareged together to produce the corresponding 
        value in the smoothed power spectrum
    """

    if start_freq:
        ind_cut = np.argmin(abs(da[dim].values-start_freq))

        spec_low = da.isel(**{dim:slice(None, ind_cut)})
        spec_high = da.isel(**{dim:slice(ind_cut, None)})
        high_smooth = spec_high.rolling(center=True, **{dim:npoints}).mean()
        spec_dani = xr.concat([spec_low, high_smooth], dim=dim)
        spec_dani = spec_dani.interpolate_na(dim=dim, method='linear')
        box_width = xr.DataArray(data=np.ones(len(spec_dani[dim])), dims=dim, coords={dim:spec_dani[dim]})
        box_width[ind_cut:] = box_width[ind_cut:]*npoints
    else:
        spec_dani = da.rolling(center=True, **{dim:npoints}).mean()
        box_width = xr.DataArray(data=np.ones(len(da[dim]))*npoints,
                dims=dim, coords={dim:da[dim]})

    return spec_dani, box_width

def spec_error(nu):
    """ Computes spectral uncertainty
    
    Parameters
    ----------
    da : array_like
        number of degrees of freedom
    
    Returns
    -------
    error_low : array_like
        lower uncertainty bound 
    error_high : array_like
        upper uncertainty bound
    """
    error_low = nu/chi2.ppf(1-0.05/2, df=nu)
    error_high = nu/chi2.ppf(0.05/2, df=nu)

    return error_low, error_high 

def wrapped_window(npts, ndata, window):
    """ Create window with the peak at the first point  
    Return an array with lenght ndata points, in which
    the window weights with npts are arranged such that the 
    the peak of the window is at the first point in the array
    and points to the left of the peak are at the 
    end of the array. 
    
    Parameters
    ----------
    npts : int
        number of window weights (needs to be odd)
    ndata : int
        number of points in the desired output array
    window : str
        boxcar or numpy window name (blackman, hanning, hamming ...)
    
    Returns
    -------
    ww : array_like
        window weights 
    """
    assert npts%2!=0, 'npts has to be odd'
    if window=='boxcar':
        w = np.ones(npts)
    else:
        w = eval('np.'+window+'(npts)')
    w = w/w.sum()
    ww = np.zeros((ndata,), dtype=float)
    mid = (npts - 1) // 2
    ww[:mid+1] = w[mid:]
    ww[-mid:] = w[:mid]
    return ww


def window_fft(npts, ndata, window):
    """ Computes the transfer function of a given window 
    
    Parameters
    ----------
    npts : int
        number of window weights (needs to be odd)
    ndata : int
        number of points in the desired output array
    window : str
        boxcar or numpy window name (blackman, hanning, hamming ...)
    
    Returns
    -------
    freq : array_like
        FFT frquencies or wavenumbers
    amp : array_like
        power spectrum of the input window
    """
    w = wrapped_window(npts, ndata, window)
    wh = np.fft.fft(w)
    ps = abs(np.fft.fftshift(wh))**2
    freq = np.fft.fftshift(np.fft.fftfreq(len(ps)))
    amp = 2*ps[freq>0]
    amp = amp/amp.max()
    freq = freq[freq>0]
    return freq, amp 


def convolve_filter(data, npts, window):
    """ Computes the convolution with a given window 
    
    Parameters
    ----------
    data : array_like
        input data to convolve
    npts : int
        number of window weights (needs to be odd)
    window : str
        boxcar or numpy window name (blackman, hanning, hamming ...)
    
    Returns
    -------
    conv : array_like
        convolution between input data and input window
    """
    
    if window=='boxcar':
        w = np.ones(npts)
    else:
        w = eval('np.'+window+'(npts)')
    w = w/w.sum()
    return np.convolve(data, w, mode="same")


def xr_convolve(da, npts, window):
    """ Convolves and `xarray.DataArray` with a given window  
    
    Parameters
    ----------
    da : `xarray.DataArray`
        input data to convolve
    npts : int
        number of window weights (needs to be odd)
    window : str
        boxcar or numpy window name (blackman, hanning, hamming ...)
    
    Returns
    -------
    da_smooth : `xarray.DataArray`
        convolution between input data and input window
    """
    tmp = []
    if 'time' in da.dims:
        for i in range(len(da.time)):         
            tmp.append(convolve_filter(da.isel(time=i).values,
                                       npts=npts, window=window))
    else:
        tmp = convolve_filter(da.values,npts=npts, window=window)
        
    da_smooth = xr.DataArray(data=np.array(tmp),
                             dims=da.dims,
                             coords=da.coords)
    return da_smooth


def filtered_spectra(da, npts, window, decimate=False):
    """ Computes the spectrum of the filterred data
     
    Computes the power spectrum of the along-track SSH after
    filtering the SSH using the window lenghts given by npts
    
    Parameters
    ----------
    da : `xarray.DataArray`
        input data (along-track SSH)
    npts : list
        list with the number of 
        window weights (each value needs to be odd)
    window : str
        boxcar or numpy window name (blackman, hanning, hamming ...)
    decimate: boolean, optional
        weather or not to subsample the data after filtering it.
        The default is False.
    
    Returns
    -------
    da_smooth : dictionary
        filtered data using a window with 
        each size in npts
        
    spec_smooth : dictionary
        spectra of the filtered data using a window with 
        each size in npts
    """
    
    dim = 'atrack_bin'
    cutoff = [1, 250, 500, 1000]

    da_smooth = {}
    spec_smooth = {}
    for count, n in enumerate(npts):
        da_convolve = xr_convolve(da, n, window=window)
        if decimate:
            tmp = da_convolve.isel(atrack_bin=slice(0,-1, decimate))
        else:
            tmp = da_convolve
        
        da_smooth['cutoff_{}'.format(cutoff[count])] = tmp
        specu = spec_xr(tmp, dim, window=True)
        specu_dani, box_width = daniell(specu,
                                       dim='freq_atrack_bin',
                                       npoints=21, start_freq=1e-4)
        spec_smooth['ssha_spec{}'.format(cutoff[count])] = specu_dani
        
    return da_smooth, spec_smooth


def add_second_axis_ssh(ax):
    """ Add a x-axis at the top of the spectra figures """
    ax2 = ax.twiny() 
    ax2.set_xscale('log')
    ax2.set_xlim(ax.axis()[0], ax.axis()[1])
    kp = np.array([0.001, 0.01, 0.1, 1.,10.,100.])
    lp = ['', '10$^2$', '10$^1$', '10$^0$', '10$^{-1}$', '10$^{-2}$']
    ax2.set_xticks(kp)
    ax2.set_xticklabels(lp)
    #plt.xlabel('Wavelength [km]', labelpad=14, fontsize=20)

def add_second_axis_window(ax):
    """ Add a x-axis at the top of the spectra figures """
    ax2 = ax.twiny() 
    ax2.set_xscale('log')
    kp = np.array([0.01, 0.1, 1.,10.,100])
    lp = ['10$^2$', '10$^1$', '10$^0$', '10$^{-1}$', '10$^{-2}$']
    ax2.set_xticks(kp)
    ax2.set_xticklabels(lp)
    #plt.xlabel('Wavelength [km]', labelpad=14, fontsize=20)
    
def add_second_axis_env(ax):
    """ Add wavelenght axis top of the spectra of env conditons"""
    ax2 = ax.twiny() 
    ax2.set_xscale('log')
    ax2.set_xlim(ax.axis()[0], ax.axis()[1])
    kp = np.array([1.,10.,100.,1000])
    lp = ['10$^0$', '10$^{-1}$', '10$^{-2}$', '10$^{-3}$']
    ax2.set_xticks(kp)
    ax2.set_xticklabels(lp)
    plt.xlabel('Wavelength [km]', labelpad=14)