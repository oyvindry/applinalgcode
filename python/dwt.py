import random
from numpy import *
import matplotlib.pyplot as plt
from sound import *
from images import *
import scipy.io as sio
import os
import os.path

def dwt_impl(x, wave_name, m = 1, bd_mode = 'symm', prefilter_mode = 'none', dims = 0, dual = False, transpose = False, data_layout = 'resolution'):
    """
    Main function for computing the DWT of a given signal. Can be used for all signals up to dimension 3.
    The dimension of the data may be one higher than the dimension of the transform, in which case the last dimension is used for 
    parallel computation.
    Note that this function computes all quantities needed from scratch in order to compute the DWT for the wavelet in question. 
    This can be time-consuming, and can be avoided by using the functions find_wav_props, find_kernel,  
    the internal DWT functions dwt1_impl_internal, dwt2_impl_internal, dwt3_impl_internal, as well as built-in persistence functions. 
    An example with minimum set of parameters is as follows:
     
    wav_props, dual_wav_props = find_wav_props(wave_name)
    save('wav_props.mat', 'wav_props', 'dual_wav_props')
    ...
    load('wav_props.mat')
    f, prefilter = find_kernel(wav_props, dual_wav_props, True)
    dwt1_impl_internal(x, f)
        
    x:         Matrix whose DWT will be computed along the first dimension(s).      
    wave_name: Name of the wavelet. Possible names are:
               'cdf97' - CDF 9/7 wavelet
               'spline53' - Spline 5/3 wavelet
               'splinex.x' - Spline wavelet with given number of vanishing moments for each filter
               'pwl0'  - Piecewise linear wavelet with 0 vanishing moments
               'pwl2'  - Piecewise linear wavelet with 2 vanishing moments
               'Haar'  - The Haar wavelet
               'dbX'   - Daubechies orthnormal wavelet with X vanishing
                         moments
               'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
                         with X vanishing moments
    m:         Number of resolutions. Default: 1.
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'none'   - Take no extra action at the boundaries
               'bd'     - Boundary wavelets
    prefilter_mode: Possible modes are:
               'none' (default)
               'filter'
               'bd_pre' - Boundary wavelets with preconditioning
    dims:      the number of dimensions to apply the DWT to. Always applied to the first dimensions. Default: max(dim(x)-1,1).
               This means that sound with many channels, and images with many colour components default to a one- and two-dimensional DWT, respectively
    dual:      Whether to apply the dual wavelet rather than the wavelet itself. Default: False
    transpose: Whether the transpose is to be taken. Default: False
    data_layout: How data should be assembled. Possible modes are:
               'resolution': Lowest resolution first (default)
               'time': Sort according to time
    """
    
    if dims == 0:
        dims = 1
        if len(shape(x)) > 1:
            dims = len(shape(x)) - 1

    wav_propsx, dual_wav_propsx = find_wav_props(wave_name, m, bd_mode, shape(x)[0])
    fx, prefilterx = find_kernel(wav_propsx, dual_wav_propsx, True, dual, transpose, prefilter_mode)
    offsets = array([[wav_propsx.offset_L, wav_propsx.offset_R]])
 
    if dims == 1:
        if transpose: # if transpose, then f will we an idwt_kernel
            idwt1_impl_internal(x, fx, m, bd_mode, prefilterx, offsets, data_layout)
        else:
            dwt1_impl_internal(x, fx, m, bd_mode, prefilterx, offsets, data_layout)
    else:
        wav_propsy, dual_wav_propsy = find_wav_props(wave_name, m, bd_mode, shape(x)[1])
        fy, prefiltery = find_kernel(wav_propsy, dual_wav_propsy, True, dual, transpose, prefilter_mode)
        offsets = concatenate( ( offsets, [[wav_propsy.offset_L, wav_propsy.offset_R]] ))
        if dims == 2:
            if transpose: # if transpose, then f will we an idwt_kernel
                idwt2_impl_internal(x, fx, fy, m, bd_mode, prefilterx, prefiltery, offsets, data_layout)
            else:
                dwt2_impl_internal(x, fx, fy, m, bd_mode, prefilterx, prefiltery, offsets, data_layout)
        else:
            wav_propsz, dual_wav_propsz = find_wav_props(wave_name, m, bd_mode, shape(x)[2])
            fz, prefilterz = find_kernel(wav_propsz, dual_wav_propsz, True, dual, transpose, prefilter_mode)
            offsets = concatenate( ( offsets, [[wav_propsz.offset_L, wav_propsz.offset_R]] ))
            if dims == 3: # if not give error message
                if transpose: # if transpose, then f will we an idwt_kernel
                    idwt3_impl_internal(x, fx, fy, fz, m, bd_mode, prefilterx, prefiltery, prefilterz, offsets, data_layout)
                else:
                    dwt3_impl_internal(x, fx, fy, fz, m, bd_mode, prefilterx, prefiltery, prefilterz, offsets, data_layout)

def idwt_impl(x, wave_name, m=1, bd_mode = 'symm', prefilter_mode = 'none', dims = 0, dual = False, transpose = False, data_layout = 'resolution'):
    """
    Main function for computing the IDWT of a given signal. Can be used for
    all signals up to dimension 3.  The dimension of the data may be one
    higher than the dimension of the transform, in which case the last
    dimension is used for parallel computation.
    
    Note that this function computes all quantities needed from scratch in
    order to compute the IDWT for the wavelet in question.  This can be
    time-consuming, and can be avoided by using the functions find_wav_props,
    find_kernel,  the internal IDWT functions idwt1_impl_internal,
    idwt2_impl_internal, idwt3_impl_internal, as well as Matlabs persistence
    functions.  An example with minimum set of parameters is as follows:
    
    [wav_props, dual_wav_props] = find_wav_props(wave_name);
    save('wav_props.mat', 'wav_props', 'dual_wav_props');
    ...
    load('wav_props.mat');
    [f, prefilter] = find_kernel(wav_props, dual_wav_props, 0);
    x = idwt1_impl_internal(x, f);
       
    x:         Matrix whose IDWT will be computed along the first dimension(s).      
    wave_name: Name of the wavelet. Possible names are:
               'cdf97' - CDF 9/7 wavelet
               'spline53' - Spline 5/3 wavelet
               'splinex.x' - Spline wavelet with given number of vanishing 
                             moments for each filter
               'pwl0'  - Piecewise linear wavelet with 0 vanishing moments
               'pwl2'  - Piecewise linear wavelet with 2 vanishing moments
               'Haar'  - The Haar wavelet
               'dbX'   - Daubechies orthnormal wavelet with X vanishing
                         moments
               'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
                         with X vanishing moments
    m:         Number of resolutions. Default: 1.
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'none'   - Take no extra action at the boundaries
               'bd'     - Boundary wavelets
    prefilter_mode: Possible modes are:
               'none' (default)
               'filter'
               'bd_pre' - Boundary wavelets with preconditioning
    dims:      the number of dimensions to apply the IDWT to. Always applied 
               to the first dimensions. Default: max(dim(x)-1,1).
               This means that sound with many channels, and images with 
               many colour components default to a one- and two-dimensional 
               IDWT, respectively
    dual:      Whether to apply the dual wavelet rather than the wavelet 
               itself. Default: 0
    transpose: Whether the transpose is to be taken. Default: 0
    data_layout: How data should be assembled. Possible modes are:
               'resolution': Lowest resolution first (default)
               'time': Sort according to time
    """
    
    if dims == 0:
        dims = 1
        if len(shape(x)) > 1:
            dims = len(shape(x)) - 1

    wav_propsx, dual_wav_propsx = find_wav_props(wave_name, m, bd_mode, shape(x)[0])
    fx, prefilterx = find_kernel(wav_propsx, dual_wav_propsx, False, dual, transpose, prefilter_mode)
    offsets = array([[wav_propsx.offset_L, wav_propsx.offset_R]])
    if dims == 1:
        if transpose: # if transpose, then f will we a dwt_kernel, 
            dwt1_impl_internal(x, fx, m, bd_mode, prefilterx, offsets, data_layout)  
        else:
            idwt1_impl_internal(x, fx, m, bd_mode, prefilterx, offsets, data_layout)
    else:
        wav_propsy, dual_wav_propsy = find_wav_props(wave_name, m, bd_mode, shape(x)[1])
        fy, prefiltery = find_kernel(wav_propsy, dual_wav_propsy, 0, dual, transpose, prefilter_mode)
        offsets = concatenate( ( offsets, [[wav_propsy.offset_L, wav_propsy.offset_R]] ))
        if dims == 2:
            if transpose: # if transpose, then f will we a dwt_kernel, 
                dwt2_impl_internal(x, fx, fy, m, bd_mode, prefilterx, prefiltery, offsets, data_layout)  
            else:
                idwt2_impl_internal(x, fx, fy, m, bd_mode, prefilterx, prefiltery, offsets, data_layout)
        else:
            wav_propsz, dual_wav_propsz = find_wav_props(wave_name, m, bd_mode, shape(x)[2])
            fz, prefilterz = find_kernel(wav_propsz, dual_wav_propsz, 0, dual, transpose, prefilter_mode)
            offsets = concatenate( ( offsets, [[wav_propsz.offset_L, wav_propsz.offset_R]] ))
            if dims == 3: # if not give error message
                if transpose: # if transpose, then f will we a dwt_kernel, 
                    dwt3_impl_internal(x, fx, fy, fz, m, bd_mode, prefilterx, prefiltery, prefilterz, offsets, data_layout)    
                else:
                    idwt3_impl_internal(x, fx, fy, fz, m, bd_mode, prefilterx, prefiltery, prefilterz, offsets, data_layout)
              


        
def dwt1_impl_internal(x, f, m = 1, bd_mode = 'symm', prefilter = 0, offsets = 0, data_layout = 'resolution'):
    """
    Compute a 1D DWT using a pre-computed kernel. The kernel may be the default library kernel obtained by calling find_kernel, 
    or may be used-defined.
    
    x:         Matrix whose DWT will be computed along the first dimension(s). 
    f:         kernel function     
    m:         Number of resolutions. Default is 1
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'none'   - Take no extra action at the boundaries
               'bd'     - Boundary wavelets
    prefilter: function which computes prefiltering. The default is no prefiltering.
    offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    data_layout: How data should be assembled. Possible modes are:
               'resolution': Lowest resolution first (default)
               'time': Sort according to time
    """
    
    if prefilter == 0:
        prefilter = lambda x, forward: x
    if len(offsets) == 0:
        offsets = array([[0,0]])
    
    prefilter(x, True)
    inds = arange(shape(x)[0])
    for res in range(m):
        xcopy = x[inds].copy()
        f(xcopy, bd_mode)
        x[inds] = xcopy
        inds = inds[offsets[0, 0]:(len(inds) - offsets[0,1]):2]
    xcopy = x[inds].copy()
    prefilter(xcopy, False)
    x[inds] = xcopy
    _reorganize_coeffs_forward(x, m, offsets, data_layout)  
    
def idwt1_impl_internal(x, f, m = 1, bd_mode = 'symm', prefilter = 0, offsets = 0, data_layout = 'resolution'):
    """
    Compute a 1D IDWT using a precomputed kernel. The kernel may be the default library kernel obtained by calling find_kernel, 
    or may be used-defined.
    
    x:         Matrix whose IDWT will be computed along the first dimension(s). 
    f:         kernel function     
    m:         Number of resolutions. Default is 1
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'none'   - Take no extra action at the boundaries
               'bd'     - Boundary wavelets
    prefilter: function which computes prefiltering. The default is no prefiltering.
    offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    data_layout: How data should be assembled. Possible modes are:
               'resolution': Lowest resolution first (default)
               'time': Sort according to time
    """
    
    if prefilter == 0:
        prefilter = lambda x, forward: x
    if len(offsets) == 0:
        offsets = array([[0,0]])
    
    resstart, resend = _reorganize_coeffs_reverse(x, m, offsets, data_layout)
    inds = arange(resstart[0,m], resend[0,m] + 1, 2**m)
    xcopy = x[inds].copy()
    prefilter(xcopy, True)
    x[inds] = xcopy
    for res in range(m-1,-1,-1):
        inds = arange(resstart[0,res], resend[0,res] + 1, 2**res)
        xcopy = x[inds].copy()
        f(xcopy, bd_mode)
        x[inds] = xcopy
    prefilter(x, False)
    
def dwt2_impl_internal(x, fx, fy, m = 1, bd_mode = 'symm', prefilterx = 0, prefiltery = 0, offsets = 0, data_layout = 'resolution'):
    """
    Compute a 2D DWT using precomputed kernels. The kernels may be the default library kernels obtained by calling find_kernel, 
    or may be used-defined.
    
    x:         Matrix whose DWT2 will be computed along the first dimensions. 
    fx, fy:    kernel functions     
    m:         Number of resolutions. Default is 1
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'none'   - Take no extra action at the boundaries
               'bd'     - Boundary wavelets
    prefilterx, prefiltery: functions which compute prefiltering. The default is no prefiltering.
    offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    data_layout: How data should be assembled. Possible modes are:
               'resolution': Lowest resolution first (default)
               'time': Sort according to time
    """
    
    if prefilterx == 0:
        prefilterx = lambda x, forward: x
    if prefiltery == 0:
        prefiltery = lambda x, forward: x
    if len(offsets) == 0:
        offsets = array([[0,0],[0,0]])
    
    # preconditioning   
    tensor2_impl(x, lambda x, bd_mode: prefilterx(x, True), lambda x, bd_mode: prefiltery(x, True), bd_mode)
         
    indsx = arange(shape(x)[0])
    indsy = arange(shape(x)[1])
    for res in range(m):
        xcopy = x[ix_(indsx, indsy)].copy()
        tensor2_impl(xcopy, fx, fy, bd_mode)
        x[ix_(indsx,indsy)] = xcopy
        indsx = indsx[offsets[0, 0]:(len(indsx) - offsets[0,1]):2]
        indsy = indsy[offsets[1, 0]:(len(indsy) - offsets[1,1]):2]
    
    # postconditioning
    xcopy = x[ix_(indsx, indsy)].copy()
    tensor2_impl(xcopy, lambda x, bd_mode: prefilterx(x, False), lambda x, bd_mode: prefiltery(x, False), bd_mode)
    x[ix_(indsx,indsy)] = xcopy
    
    _reorganize_coeffs2_forward(x, m, offsets, data_layout)
    
def idwt2_impl_internal(x, fx, fy, m = 1, bd_mode = 'symm', prefilterx = 0, prefiltery = 0, offsets = 0, data_layout = 'resolution'):
    """
    Compute a 2D IDWT using pre-computed kernels. The kernels may be the default library kernels obtained by calling find_kernel, 
    or may be used-defined.
    
    x:         Matrix whose IDWT2 will be computed along the first dimensions. 
    fx, fy:    kernel functions     
    m:         Number of resolutions. Default is 1
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'none'   - Take no extra action at the boundaries
               'bd'     - Boundary wavelets
    prefilterx, prefiltery: functions which compute prefiltering. The default is no prefiltering.
    offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    data_layout: How data should be assembled. Possible modes are:
               'resolution': Lowest resolution first (default)
               'time': Sort according to time
    """
    
    if prefilterx == 0:
        prefilterx = lambda x, forward: x
    if prefiltery == 0:
        prefiltery = lambda x, forward: x
    if len(offsets) == 0:
        offsets = array([[0,0],[0,0]])
    
    resstart, resend = _reorganize_coeffs2_reverse(x, m, offsets, data_layout)
        
    # postconditioning
    indsx = arange(resstart[0,m], resend[0,m] + 1, 2**m)
    indsy = arange(resstart[1,m], resend[1,m] + 1, 2**m)
    xcopy = x[ix_(indsx, indsy)].copy()
    tensor2_impl(xcopy, lambda x, bd_mode: prefilterx(x, True), lambda x, bd_mode: prefiltery(x, True), bd_mode)
    x[ix_(indsx, indsy)] = xcopy

    for res in range(m - 1, -1, -1):
        indsx = arange(resstart[0, res], resend[0, res] + 1, 2**res)
        indsy = arange(resstart[1, res], resend[1, res] + 1, 2**res)
        xcopy = x[ix_(indsx, indsy)].copy()
        tensor2_impl(xcopy, fx, fy, bd_mode)
        x[ix_(indsx, indsy)] = xcopy
    
    # preconditioning
    indsx = arange(resstart[0, 0], resend[0, 0] + 1)
    indsy = arange(resstart[1, 0], resend[1, 0] + 1)
    xcopy = x[ix_(indsx, indsy)].copy()
    tensor2_impl(xcopy, lambda x, bd_mode: prefilterx(x, False), lambda x, bd_mode: prefiltery(x, False), bd_mode)
    x[ix_(indsx, indsy)] = xcopy

def dwt3_impl_internal(x, fx, fy, fz, m = 1, bd_mode = 'symm', prefilterx = 0, prefiltery = 0, prefilterz = 0, offsets = 0, data_layout = 'resolution'):
    """
    Compute a 3D DWT using precomputed kernels. The kernels may be the default library kernels obtained by calling find_kernel, 
    or may be used-defined.
        % x:         Matrix whose DWT3 will be computed along the first dimensions. 
    fx, fy, fz: kernel functions     
    m:         Number of resolutions. Default is 1
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'none'   - Take no extra action at the boundaries
               'bd'     - Boundary wavelets
    prefilterx, prefiltery, prefilterz: functions which compute prefiltering. The default is no prefiltering.
    offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    data_layout: How data should be assembled. Possible modes are:
               'resolution': Lowest resolution first (default)
               'time': Sort according to time
    """
    
    if prefilterx == 0:
        prefilterx = lambda x, forward: x
    if prefiltery == 0:
        prefiltery = lambda x, forward: x
    if prefilterz == 0:
        prefilterz = lambda x, forward: x
    if len(offsets) == 0:
        offsets = array([[0,0],[0,0],[0,0]]) 

    # preconditioning   
    tensor3_impl(x, lambda x, bd_mode: prefilterx(x, True), lambda x, bd_mode: prefiltery(x, True), lambda x, bd_mode: prefilterz(x, True), bd_mode)
    
    indsx = arange(shape(x)[0])
    indsy = arange(shape(x)[1])
    indsz = arange(shape(x)[2])
    for res in range(m):
        xcopy = x[ix_(indsx, indsy, indsz)].copy()
        tensor3_impl(xcopy, fx, fy, fz, bd_mode)
        x[ix_(indsx, indsy, indsz)] = xcopy
        indsx = indsx[offsets[0, 0]:(len(indsx) - offsets[0,1]):2]
        indsy = indsy[offsets[1, 0]:(len(indsy) - offsets[1,1]):2]
        indsz = indsz[offsets[2, 0]:(len(indsz) - offsets[2,1]):2]
    
    # postconditioning
    xcopy = x[ix_(indsx, indsy, indsz)].copy()
    tensor3_impl(xcopy, lambda x, bd_mode: prefilterx(x, False), lambda x, bd_mode: prefiltery(x, False), lambda x, bd_mode: prefilterz(x, False), bd_mode)
    x[ix_(indsx, indsy, indsz)] = xcopy
    
    _reorganize_coeffs3_forward(x, m, offsets, data_layout)
        
        
    
    
    
    
    
def idwt3_impl_internal(x, fx, fy, fz, m = 1, bd_mode = 'symm', prefilterx = 0, prefiltery = 0, prefilterz = 0, offsets = 0, data_layout = 'resolution'):
    """
    Compute a 3D IDWT using precomputed kernels. The kernels may be the default library kernels obtained by calling find_kernel, 
    or may be used-defined.
    
    x:         Matrix whose IDWT3 will be computed along the first dimensions. 
    fx, fy, fz: kernel functions     
    m:         Number of resolutions. Default is 1
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'none'   - Take no extra action at the boundaries
               'bd'     - Boundary wavelets
    prefilterx, prefiltery, prefilterz: functions which compute prefiltering. The default is no prefiltering.
    offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    data_layout: How data should be assembled. Possible modes are:
               'resolution': Lowest resolution first (default)
               'time': Sort according to time
    """
    
    if prefilterx == 0:
        prefilterx = lambda x, forward: x
    if prefiltery == 0:
        prefiltery = lambda x, forward: x
    if prefilterz == 0:
        prefilterz = lambda x, forward: x
    if len(offsets) == 0:
        offsets = array([[0,0],[0,0],[0,0]])
    
    resstart, resend = _reorganize_coeffs3_reverse(x, m, offsets, data_layout)
    
    # postconditioning
    indsx = arange(resstart[0,m], resend[0,m] + 1, 2**m)
    indsy = arange(resstart[1,m], resend[1,m] + 1, 2**m)
    indsz = arange(resstart[2,m], resend[2,m] + 1, 2**m)
    xcopy = x[ix_(indsx, indsy, indsz)].copy()
    tensor3_impl(xcopy, lambda x, bd_mode: prefilterx(x, True), lambda x, bd_mode: prefiltery(x, True), lambda x, bd_mode: prefilterz(x, True), bd_mode)
    x[ix_(indsx, indsy, indsz)] = xcopy
    
    for res in range(m - 1, -1, -1):
        indsx = arange(resstart[0, res], resend[0, res] + 1, 2**res)
        indsy = arange(resstart[1, res], resend[1, res] + 1, 2**res)
        indsz = arange(resstart[2, res], resend[2, res] + 1, 2**res)
        xcopy = x[ix_(indsx, indsy, indsz)].copy()
        tensor3_impl(xcopy, fx, fy, fz, bd_mode)
        x[ix_(indsx, indsy, indsz)] = xcopy
    
    # preconditioning
    indsx = arange(resstart[0, 0], resend[0, 0] + 1)
    indsy = arange(resstart[1, 0], resend[1, 0] + 1)
    indsz = arange(resstart[2, 0], resend[2, 0] + 1)
    xcopy = x[ix_(indsx, indsy, indsz)].copy()
    tensor3_impl(xcopy, lambda x, bd_mode: prefilterx(x, False), lambda x, bd_mode: prefiltery(x, False), lambda x, bd_mode: prefilterz(x, False), bd_mode)                      
    x[ix_(indsx, indsy, indsz)] = xcopy
    
        
    
    

    
    
    
    
    
    
    
    
    
    
    
class WaveletProps(object):
    def __init__(self, wave_name, m, length_signal, offset_L, offset_R):
        self.wave_name = wave_name
        self.m = m
        self.length_signal = length_signal
        self.offset_L = 0
        self.offset_R = 0
        
    
def find_wav_props(wave_name, m = 1, bd_mode = 'symm', length_signal = 0):
    """
    Computes the properties of a wavelet with the given name. What properties are computed depend on the bd_mode parameter, m, and length_signal.
    
    wave_name: Name of the wavelet. Possible names are:
               'cdf97' - CDF 9/7 wavelet
               'spline53' - Spline 5/3 wavelet
               'splinex.x' - Spline wavelet with given number of vanishing moments for each filter
               'pwl0'  - Piecewise linear wavelet with 0 vanishing moments
               'pwl2'  - Piecewise linear wavelet with 2 vanishing moments
               'Haar'  - The Haar wavelet
               'dbX'   - Daubechies orthonormal wavelet with X vanishing
                         moments
               'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
                         with X vanishing moments
    m:         Number of resolutions. Default: 1
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'none'   - Take no extra action at the boundaries
               'bd'     - Boundary wavelets
    length_signal: Length of the input signal. Default: 0.
    """
    
    wav_props = WaveletProps(wave_name, m, length_signal, 0, 0)
    dual_wav_props = WaveletProps(wave_name, m, length_signal, 0, 0)
    
    if wave_name.lower()[:2] == 'db':
        N = int(wave_name[2:])
        if N > 1:
            _wav_props_ortho(N, wav_props, dual_wav_props, bd_mode)
    elif wave_name.lower()[:3] == 'sym':
        N = int(wave_name[3:]);
        if N > 1:
            _wav_props_ortho(N, wav_props, dual_wav_props, bd_mode, 1)
    elif wave_name.lower() == 'pwl0':
        _wav_props_pwl0(wav_props, dual_wav_props, bd_mode)
    elif wave_name.lower() == 'pwl2':
        _wav_props_pwl2(wav_props, dual_wav_props, bd_mode)
    elif wave_name.lower() == 'spline53':
        _wav_props_53(wav_props, dual_wav_props, bd_mode)
    elif wave_name.lower() == 'cdf97':
        _wav_props_97(wav_props, dual_wav_props, bd_mode)
    elif wave_name.lower()[:6] =='spline':
        N = int(wave_name[6])
        Ntilde = int(wave_name[8])
        dual_wav_props.g0, dual_wav_props.g1, wav_props.g0, wav_props.g1 = _compute_spline_filters(N, Ntilde)
        if bd_mode.lower() == 'bd':
            WL, WLtilde, WR, WRtilde = _wav_props_biortho_bd(N, Ntilde, wav_props, dual_wav_props)
            # swap filters if offset is odd
            if mod(wav_props.offset_L,2) == 1:
                g0temp = wav_props.g0; wav_props.g0 = wav_props.g1; wav_props.g1 = g0temp
                g0temp = dual_wav_props.g0; dual_wav_props.g0 = dual_wav_props.g1; dual_wav_props.g1 = g0temp
            find_AL_AR(WL, WR, wav_props)
            find_AL_AR(WLtilde, WRtilde, dual_wav_props)
    return wav_props, dual_wav_props


            
def _wav_props_ortho(N, wav_props, dual_wav_props, bd_mode, type = 0): 
    """
    N:    Number of vanishing moments
    type: The type of orthonormal wavelet.
          0: Daubechies wavelets with minimum phase (default)  
          1: Symmlets - wavelets with close to linear phase (almost symmetric)
    
    """
        
    h0 = 0; h1 = 0
    if type == 0:
        [h0, h1, wav_props.g0, wav_props.g1] = _h0h1_compute_ortho(N)
    elif type == 1:
        [h0, h1, wav_props.g0, wav_props.g1] = _h0h1_compute_sym(N)
    wav_props.lambdas, wav_props.alpha, wav_props.beta, wav_props.last_even = _lifting_fact_ortho(h0, h1)
    
    # TODO: add boundary handling
    
    dual_wav_props.lambdas = -fliplr(wav_props.lambdas)
    dual_wav_props.alpha = 1/wav_props.alpha
    dual_wav_props.beta = 1/wav_props.beta
    dual_wav_props.last_even = not wav_props.last_even
    
    # TODO: add boundary handling
    
    
    
def _set_wav_props(wav_props, dual_wav_props, bd_mode): # TODO: Add WL, WLtilde, WR, WRtilde
    # TODO: add boundary handling
    
    dual_wav_props.lambdas = -wav_props.lambdas 
    dual_wav_props.alpha = 1/wav_props.alpha 
    dual_wav_props.beta = 1/wav_props.beta
    dual_wav_props.last_even = not wav_props.last_even
    
    # TODO: add boundary handling

            
def _wav_props_pwl0(wav_props, dual_wav_props, bd_mode):
    wav_props.last_even = False
    wav_props.lambdas = array([0.5])
    wav_props.alpha = sqrt(2)
    wav_props.beta = sqrt(2)
    dual_wav_props.g0 = [sqrt(2)]
    dual_wav_props.g1 = array([-0.5, 1., -0.5])*sqrt(2)
    wav_props.g0 = [1/2., 1, 1/2.]/sqrt(2)
    wav_props.g1 = [1/sqrt(2)]
    
    # TODO: add boundary handling
    
    _set_wav_props(wav_props, dual_wav_props, bd_mode); # TODO: Add WL, WLtilde, WR, WRtilde

def _wav_props_pwl2(wav_props, dual_wav_props, bd_mode):
    wav_props.last_even = False
    wav_props.lambdas = array([0.25, 0.5])
    wav_props.alpha = sqrt(2)
    wav_props.beta = sqrt(2)
    dual_wav_props.g0 = array([-1/8., 1/4., 3/4., 1/4., -1/8.])*sqrt(2)
    dual_wav_props.g1 = array([-1/2., 1, -1/2.])*sqrt(2)
    wav_props.g0 = array([1/2., 1, 1/2.])/sqrt(2)
    wav_props.g1 = array([-1/8., -1/4., 3/4., -1/4., -1/8.])/sqrt(2)
    
    # TODO: add boundary handling
    
    _set_wav_props(wav_props, dual_wav_props, bd_mode); # TODO: Add WL, WLtilde, WR, WRtilde

def _wav_props_53(wav_props, dual_wav_props, bd_mode):
    wav_props.last_even = False;
    wav_props.lambdas = array([-1, 0.125])
    wav_props.alpha = 2
    wav_props.beta = 0.5
    dual_wav_props.g0 = array([-1/4, 1/2, 3/2, 1/2, -1/4])
    dual_wav_props.g1 = array([-1/4, 1/2, -1/4])
    wav_props.g0 = array([1/4, 1/2, 1/4])
    wav_props.g1 = array([-1/4, -1/2, 3/2, -1/2, -1/4])
    
    # TODO: add boundary handling
    
    _set_wav_props(wav_props, dual_wav_props, bd_mode); # TODO: Add WL, WLtilde, WR, WRtilde

def _wav_props_97(wav_props, dual_wav_props, bd_mode):
    wav_props.last_even = False;
    wav_props.lambdas, wav_props.alpha, wav_props.beta, dual_wav_props.g0, dual_wav_props.g1, wav_props.g0, wav_props.g1 = _lifting_fact_97()
    
    # TODO: add boundary handling
    
    _set_wav_props(wav_props, dual_wav_props, bd_mode); # TODO: Add WL, WLtilde, WR, WRtilde
            



def _lifting_fact_97():
    h0, h1, g0, g1 = _h0h1_compute_97() # Should have 9 and 7 filter coefficients.
    h00, h01 = h0[0:9:2], h0[1:9:2]
    h10, h11 = h1[0:7:2], h1[1:7:2]
        
    lambdas=array(zeros(4))
        
    lambdas[0] = -h00[0]/h10[0]
    h00[0:5] = h00[0:5]+convolve(h10[0:4],[lambdas[0],lambdas[0]])
    h01[0:4] = h01[0:4]+convolve(h11[0:3],[lambdas[0],lambdas[0]])  
        
    lambdas[1] = -h10[0]/h00[1]
    h10[0:4] = h10[0:4]+convolve(h00[1:4],[lambdas[1],lambdas[1]])
    h11[0:3] = h11[0:3]+convolve(h01[1:3],[lambdas[1],lambdas[1]]) 
        
    lambdas[2] = -h00[1]/h10[1]
    h00[1:4] = h00[1:4]+convolve(h10[1:3],[lambdas[2],lambdas[2]])
    h01[1:3] = h01[1:3]+convolve(h11[1:2],[lambdas[2],lambdas[2]])  
        
    lambdas[3] = -h10[1]/h00[2]
    h10[0:4] = h10[0:4]+convolve(h00[1:4],[lambdas[3],lambdas[3]])
    h11[0:3] = h11[0:3]+convolve(h01[1:3],[lambdas[3],lambdas[3]]) 
        
    alpha, beta = h00[2], h11[1] 
    return lambdas, alpha, beta,  h0, h1, g0, g1
    
def _h0h1_compute_97():
    QN = _compute_QN(4)
        
    rts = roots(QN)
    rts1 = rts[nonzero(abs(imag(rts))>0.001)] # imaginary roots
    rts2 = rts[nonzero(abs(imag(rts))<0.001)] # real roots
        
    h0=array([1])
    for rt in rts1:
        h0 = convolve(h0, [-rt,1])
    for k in range(2):
        h0 = convolve(h0,[1/4.,1/2.,1/4.])
    h0=h0*QN[0]
        
    g0=array([1])
    for rt in rts2:
        g0 = convolve(g0, [-rt,1])
    for k in range(2):
        g0 = convolve(g0,[1/4.,1/2.,1/4.])
        
    g0, h0 = real(g0), real(h0)
    x = sqrt(2)/abs(sum(h0))
    g0, h0 = g0/x, h0*x
    N= shape(g0)[0]
    h1=g0*(-1)**(arange(-(N-1)/2,(N+1)/2))
    N= shape(h0)[0]
    g1=h0*(-1)**(arange(-(N-1)/2,(N+1)/2))
    #print(h0, h1, g0, g1)
    return h0, h1, g0, g1

def _compute_QN(N):
    """
    Compute the coefficients of the polynomial Q^(N)((1-cos(w))/2).
    """
    QN=zeros(N)
    for k in range(N):
        QN[k] = 2*math.factorial(N+k-1)/(math.factorial(k)*math.factorial(N-1))
    vals = array([QN[0]])
    start = array([1.0])
    for k in range(1, N):
        start = convolve(start,[-1/4.0,1/2.0,-1/4.0])
        vals = hstack([0,vals])
        vals = hstack([vals,0])
        vals = vals + QN[k]*start
    return vals
        
def _compute_spline_filters(N, Ntilde):
    Navg=int((N+Ntilde)/2)
    vals = _compute_QN(Navg)
    
    h0=array([1])
    for k in range(int(N/2)):
        h0 = convolve(h0,[1/4.,1/2.,1/4.])
    h0 = convolve(h0, vals)
    
    g0=array([1])
    for k in range(int(Ntilde/2)):
        g0 = convolve(g0,[1/4.,1/2.,1/4.])
    
    h1=g0*(-1)**(array(range(-int((len(g0)-1)/2),int((len(g0)+1)/2))))
    g1=h0*(-1)**(array(range(-int((len(h0)-1)/2),int((len(h0)+1)/2))))
    return h0, h1, g0, g1

        
def _h0h1_compute_ortho(N):
    vals = _compute_QN(N)
    rts=roots(vals)
    rts1=rts[nonzero(abs(rts)>1)]
    g0=array([1])
    for rt in rts1:
        g0=convolve(g0,[-rt,1])
        
    K=sqrt(vals[0]*(-1)**(len(rts1))/prod(rts1))
    g0=K*g0
    for k in range(N):
        g0=convolve(g0,[1/2.,1/2.])
        
    g0=real(g0)
    h0=g0[::-1]
    g1=g0[::-1]*(-1)**(array(range(len(g0))))
    h1=g1[::-1]
    return h0, h1, g0, g1

def _lifting_fact_ortho(h0, h1):
    """
    Assume that len(h1)==len(h0), and that h0 and h1 are even length and as symmetric as possible, with h0 with a minimum possible overweight of filter coefficients to the left, h1 to the right
    This function computes lifting steps l1, l2,...,ln, and constants alpha, beta so that ln ... l2 l1 H =  diag(alpha,beta), and stores these in files.
    This gives the following recipes for 
        Computing H: first multiply with diag(alpha,beta), then the inverses of the lifting steps in reverse order 
        Computing G: apply the lifting steps in the order given, finally multiply with diag(1/alpha,1/beta)
    ln is always odd, so that l1 is odd if and only if n is odd.
    All even lifting steps have only filter coefficients 0,1. All odd lifting steps have only filter coefficients -1,0
    """
        
    stepnr=0
    start1, end1, len1, start2, end2, len2 = 0, int(len(h0)/2)-1, int(len(h0)/2),  0, int(len(h1)/2)-1, int(len(h1)/2)
    lambdas=zeros((len1+1,2))
    if mod(len1,2)==0: # Start with an even step
        h00, h01 = h0[0:len(h0):2], h0[1:len(h0):2]
        h10, h11 = h1[0:len(h1):2], h1[1:len(h1):2]
      
        lambda1=-h00[0]/h10[0]
        h00=h00+lambda1*h10 
        h01=h01+lambda1*h11
        start1, end1, len1 = 1, len1-1, len1-1
        lambdas[stepnr,:] = [lambda1,0]
    else: # Start with an odd step
        h00, h01 = h0[1:len(h0):2], h0[0:len(h0):2]
        h10, h11 = h1[1:len(h1):2], h1[0:len(h1):2]
        
        lambda1=-h10[end1]/h00[end1] 
        h10=h10+lambda1*h00 
        h11=h11+lambda1*h01
        start2, end2, len2 = 0, len2 - 2, len2-1
        lambdas[stepnr,:] = [0,lambda1]
      
     
    stepnr=stepnr+1
    
    # print([h00 h01; h10 h11], convolve(h00,h11)-convolve(h10,h01))
    while len2>0: # Stop when the second element in the first column is zero
        if len1>len2: # Reduce the degree in the first row. 
            lambda1=-h00[start1]/h10[start2]
            lambda2=-h00[end1]/h10[end2]
            h00[start1:(end1+1)] = h00[start1:(end1+1)]+convolve(h10[start2:(end2+1)],[lambda1,lambda2])
            h01[start1:(end1+1)] = h01[start1:(end1+1)]+convolve(h11[start2:(end2+1)],[lambda1,lambda2])
            start1, end1, len1 = start1+1, end1-1, len1-2
        else: # reduce the degree in the second row. 
            lambda1=-h10[start2]/h00[start1]
            lambda2=-h10[end2]/h00[end1]
            h10[start2:(end2+1)] = h10[start2:(end2+1)]+convolve(h00[start1:(end1+1)],[lambda1,lambda2])
            h11[start2:(end2+1)] = h11[start2:(end2+1)]+convolve(h01[start1:(end1+1)],[lambda1,lambda2])
            start2, end2, len2 = start2+1, end2-1, len2-2
        lambdas[stepnr,:]=[lambda1,lambda2]
        stepnr=stepnr+1
        
    # print([h00 h01; h10 h11], convolve(h00,h11)-convolve(h10,h01))
      
    # Add the final lifting, and compute alpha,beta
    alpha=sum(h00)
    beta=sum(h11)
    lastlift=-sum(h01)/beta
    if mod(len(h0)/2,2)==0:
        lambdas[stepnr,:] = [0,lastlift]
    else:
        lambdas[stepnr,:] = [lastlift,0]
    # [h00 h01; h10 h11]
    return lambdas, alpha, beta, True
    
    
    
    
        
    
# Kernel functions    
            
def find_kernel(wav_props, dual_wav_props, forward, dual = False, transpose = False, prefilter_mode = 'none'):
    """ 
    Function which returns the default kernel of the library for use with the wavelet with properties encapsulated by the given parameters.
    The kernel can be passed to the the internal functions (i)dwt1_impl_internal, (i)dwt2_impl_internal, (i)dwt3_impl_internal. 
    User-defined kernels can also be passed to these internal functions: They simply have to take the x and bd_mode parameters, and return the 
    transformed vectors.
    
    wav_props: Object which encapsulates the wavelet
    dual_wav_props: Object which encapsulates the dual wavelet
    forward: Whether the forward transform should be used
    dual: (optional). Whether the dual wavelet should be applied. Default is 0
    transpose: (optional). Default is 0
    prefilter_mode: (optional). Default is 'none'
    """
    
    prefilter = lambda x, forward: x
    if transpose:
        forward = not forward
        dual = not dual
    if dual:
        temp = wav_props; wav_props = dual_wav_props; dual_wav_props = temp
    
    if wav_props.wave_name.lower() == 'haar':
        if forward:
            f = _dwt_kernel_haar
        else:
            f = _idwt_kernel_haar
    else:
       if forward:
           f, prefilter = _find_kernel_dwt_general(wav_props, dual_wav_props, prefilter_mode)
       else:
           f, prefilter = _find_kernel_idwt_general(wav_props, dual_wav_props, prefilter_mode)
           
    return f, prefilter      
            

def _find_kernel_dwt_general(wav_props, dual_wav_props, prefilter_mode):
    prefilter = lambda x, forward: x
    
    if wav_props.wave_name.lower() == 'spline53' or wav_props.wave_name.lower() == 'cdf97' or wav_props.wave_name.lower() ==  'pwl0' or wav_props.wave_name.lower() ==  'pwl2':
        f = lambda x, bd_mode: _dwt_kernel_biortho(x, bd_mode, dual_wav_props)
        # TODO: Add boundary handling
    elif wav_props.wave_name.lower()[:2] == 'db':
        N = int(wav_props.wave_name[2:])
        if N == 1:
            f = _dwt_kernel_haar
        else:
            f = lambda x, bd_mode: _dwt_kernel_ortho(x, bd_mode, dual_wav_props)
            # TODO: Add boundary handling
    elif wav_props.wave_name.lower()[:3] == 'sym':
        N = int(wav_props.wave_name[3:])
        if N == 1:
            f = _dwt_kernel_haar
        else:
            f = lambda x, bd_mode: _dwt_kernel_ortho(x, bd_mode, dual_wav_props)
            # TODO: Add boundary handling
    elif wav_props.wave_name.lower()[:6] == 'spline':
        N = int(wav_props.wave_name[6])
        Ntilde = int(wav_props.wave_name[8])
        f = lambda x, bd_mode: _dwt_kernel_filters(x, bd_mode, dual_wav_props)
        # TODO: Add boundary handling
        
    return f, prefilter

def _find_kernel_idwt_general(wav_props, dual_wav_props, prefilter_mode):
    prefilter = lambda x, forward: x
    
    if wav_props.wave_name.lower() == 'spline53' or wav_props.wave_name.lower() == 'cdf97' or wav_props.wave_name.lower() ==  'pwl0' or wav_props.wave_name.lower() ==  'pwl2':
        f = lambda x, bd_mode: _idwt_kernel_biortho(x, bd_mode, wav_props)
        # TODO: Add boundary handling
    elif wav_props.wave_name.lower()[:2] == 'db':
        N = int(wav_props.wave_name[2:])
        if N == 1:
            f = _idwt_kernel_haar;
        else:
            f = lambda x, bd_mode: _idwt_kernel_ortho(x, bd_mode, wav_props)
            # TODO: Add boundary handling
    elif wav_props.wave_name.lower()[:3] == 'sym':
        N = int(wav_props.wave_name[3:])
        if N == 1:
            f = _idwt_kernel_haar
        else:
            f = lambda x, bd_mode: _idwt_kernel_ortho(x, bd_mode, wav_props)
            # TODO: Add boundary handling
    elif wav_props.wave_name.lower()[:6] == 'spline':
        N = int(wav_props.wave_name[6])
        Ntilde = int(wav_props.wave_name[8])
        f = lambda x, bd_mode: _idwt_kernel_filters(x, bd_mode, wav_props)
        # TODO: Add boundary handling

    return f, prefilter

# TODO Add boundary handling at the beginning and end 

def _dwt_kernel_filters(x, bd_mode, dual_wav_props):
    x0 = x.copy()
    x1 = x.copy()
    filter_impl(dual_wav_props.g0, x0, bd_mode)
    filter_impl(dual_wav_props.g1, x1, bd_mode)
    x[::2] = x0[::2]
    x[1::2] = x1[1::2]
# End _dwt_kernel_filters

# TODO Add boundary handling at the beginning and end 

def _idwt_kernel_filters(x, bd_mode, wav_props):
    x0 = x.copy(); x0[1::2] = 0
    x1 = x.copy(); x1[::2] = 0
    filter_impl(wav_props.g0, x0, bd_mode)
    filter_impl(wav_props.g1, x1, bd_mode)
    x[:] = x0 + x1
# End _idwt_kernel_filters
    
# The Haar wavelet

def _dwt_kernel_haar(x, bd_mode):
    x /= sqrt(2)
    if mod(len(x), 2)==1:
        a, b = x[0] + x[1] - x[-1], x[0] - x[1] - x[-1]
        x[0], x[1] = a, b 
        x[-1] *= 2
    else:
        a, b = x[0] + x[1], x[0] - x[1] 
        x[0], x[1] = a, b 
    for k in range(2,len(x) - 1,2):
        a, b = x[k] + x[k+1], x[k] - x[k+1]  
        x[k], x[k+1] = a, b
# End _dwt_kernel_haar
         
def _idwt_kernel_haar(x, bd_mode):
    x /= sqrt(2)
    if mod(len(x), 2)==1:
        a, b = x[0] + x[1]  + x[-1], x[0] - x[1]
        x[0], x[1] = a, b
        for k in range(2,len(x) - 2, 2):
            a, b = x[k] + x[k+1], x[k] - x[k+1] 
            x[k], x[k+1] = a, b
    else:    
        for k in range(0,len(x) - 1, 2):
            a, b = x[k] + x[k+1], x[k] - x[k+1] 
            x[k], x[k+1] = a, b
# End _idwt_kernel_haar

def _dwt_kernel_biortho(x, bd_mode, dual_wav_props):
    # TODO: Add boundary handling
    x[0::2] /= dual_wav_props.alpha
    x[1::2] /= dual_wav_props.beta
    iseven = not dual_wav_props.last_even
    for stepnr in range(dual_wav_props.lambdas.shape[0] - 1, -1, -1):
        if iseven:
            lifting_even_symm(  dual_wav_props.lambdas[stepnr], x, bd_mode)
        else:
            lifting_odd_symm(  dual_wav_props.lambdas[stepnr], x, bd_mode)
        iseven = not iseven
    # TODO: Add boundary handling
# End _dwt_kernel_biortho

def _idwt_kernel_biortho(x, bd_mode, wav_props):
    # TODO: Add boundary handling
    iseven = ( mod(wav_props.lambdas.shape[0], 2) == wav_props.last_even )
    for stepnr in range(wav_props.lambdas.shape[0]):
        if iseven:
            lifting_even_symm( wav_props.lambdas[stepnr], x, bd_mode)
        else:    
            lifting_odd_symm(  wav_props.lambdas[stepnr], x, bd_mode)
        iseven = not iseven
    x[0::2] /= wav_props.alpha
    x[1::2] /= wav_props.beta
    # TODO: Add boundary handling
# End _idwt_kernel_biortho

def _dwt_kernel_ortho(x, bd_mode, dual_wav_props):
    # TODO: Add boundary handling
    x[0::2] /= dual_wav_props.alpha
    x[1::2] /= dual_wav_props.beta
    iseven = not dual_wav_props.last_even
    for stepnr in range(dual_wav_props.lambdas.shape[0] - 1, -1, -1):
        if iseven:
            lifting_even( dual_wav_props.lambdas[stepnr, 1], \
                          dual_wav_props.lambdas[stepnr, 0], \
                          x, bd_mode)
        else:    
            lifting_odd(  dual_wav_props.lambdas[stepnr, 1], \
                          dual_wav_props.lambdas[stepnr, 0], \
                          x, bd_mode)
        iseven = not iseven
    # TODO: Add boundary handling
# End _dwt_kernel_ortho

def _idwt_kernel_ortho(x, bd_mode, wav_props):
    # TODO: Add boundary handling
    iseven = ( mod(wav_props.lambdas.shape[0], 2) == wav_props.last_even )
    for stepnr in range(wav_props.lambdas.shape[0]):
        if iseven:
            lifting_even( wav_props.lambdas[stepnr, 0], \
                          wav_props.lambdas[stepnr, 1], \
                          x, bd_mode)
        else:    
            lifting_odd(  wav_props.lambdas[stepnr, 0], \
                          wav_props.lambdas[stepnr, 1], \
                          x, bd_mode)
        iseven = not iseven
    x[0::2] /= wav_props.alpha
    x[1::2] /= wav_props.beta
    # TODO: Add boundary handling
# End _idwt_kernel_ortho

# Lifting steps
               
def lifting_even(lmbda1, lmbda2, x, bd_mode):
    if mod(len(x), 2)!=0:
        raise AssertionError()
    x[0] += lmbda1*x[1] + lmbda2*x[-1]
    x[2:-1:2] += lmbda1*x[3::2] + lmbda2*x[1:-2:2]
# End lifting_even
                
def lifting_odd(lmbda1, lmbda2, x, bd_mode):
    if mod(len(x), 2)!=0:
        raise AssertionError()
    x[1:-2:2] += lmbda1*x[2:-1:2] + lmbda2*x[0:-3:2]
    x[-1] += lmbda1*x[0] + lmbda2*x[-2]
# End lifting_odd
    
def lifting_even_symm(lmbda, x, bd_mode):
    if (bd_mode.lower() == 'per') and mod(len(x), 2)!=0:
        raise AssertionError()
    if bd_mode.lower() == 'symm':
        x[0] += 2*lmbda*x[1] # With symmetric extension
    else:
        x[0] += lmbda*(x[1]+x[-1])
    x[2:-1:2] += lmbda*(x[1:-2:2] + x[3::2])
    if mod(len(x), 2)==1 and bd_mode.lower() == 'symm':
        x[-1] += 2*lmbda*x[-2] # With symmetric extension
# End lifting_even_symm
      
def lifting_odd_symm(lmbda, x, bd_mode):
    if (bd_mode.lower() == 'per') and mod(len(x), 2) != 0:
        raise AssertionError()
    x[1:-1:2] += lmbda*(x[0:-2:2] + x[2::2])
    if mod(len(x), 2) == 0:
        if bd_mode.lower() == 'symm':
            x[-1] += 2*lmbda*x[-2] # With symmetric extension
        else:
            x[-1] += lmbda*(x[0]+x[-2])
# End lifting_odd_symm

def _reorganize_coeffs_forward(x, m, offsets, data_layout):
    if data_layout.lower() == 'resolution':
        N = shape(x)[0]
        y = zeros_like(x)
        inds = arange(N)
        endy = N
        for res in range(1, m + 1):
            xindices = concatenate( ( inds[:offsets[0,0]], inds[(offsets[0, 0] + 1):(len(inds) - offsets[0, 1]):2], inds[(len(inds) - offsets[0, 1]):] ) ) # psi-indices
            y[(endy-len(xindices)):endy] = x[xindices]
            endy = endy-len(xindices)
            inds = inds[ offsets[0,0]:(len(inds) - offsets[0,1]):2 ]
        y[:endy] = x[inds]
        x[:] = y[:]

def _reorganize_coeffs2_forward(sig_in, m, offsets, data_layout):
    if data_layout.lower() == 'resolution':
        sig_out = zeros_like(sig_in)
        indsx = arange(shape(sig_in)[0])
        indsy = arange(shape(sig_in)[1])
        endx = shape(sig_in)[0]; endy = shape(sig_in)[1]
        for res in range(1, m + 1):
            szx = len(indsx); szy = len(indsy)
            psiinds_x = concatenate( ( indsx[:offsets[0,0]], indsx[(offsets[0,0] + 1):(szx - offsets[0,1]):2], indsx[(szx - offsets[0,1]):] ) ) # psi-indices
            psiinds_y = concatenate( ( indsy[:offsets[1,0]], indsy[(offsets[1,0] + 1):(szy - offsets[1,1]):2], indsy[(szy - offsets[1,1]):] ) )          
            phiinds_x = indsx[offsets[0,0]:(szx - offsets[0,1]):2]
            phiinds_y = indsy[offsets[1,0]:(szy - offsets[1,1]):2]
            
            sig_out[ (endx-len(psiinds_x)):endx, :(endy-len(psiinds_y))] =     sig_in[ ix_(psiinds_x, phiinds_y) ]
            sig_out[ :(endx-len(psiinds_x)), (endy-len(psiinds_y)):endy] =     sig_in[ ix_(phiinds_x, psiinds_y) ]
            sig_out[ (endx-len(psiinds_x)):endx, (endy-len(psiinds_y)):endy] = sig_in[ ix_(psiinds_x, psiinds_y) ]
            
            endx -= len(psiinds_x); endy -= len(psiinds_y)
            indsx = indsx[ offsets[0,0]:(szx - offsets[0,1]):2 ] 
            indsy = indsy[ offsets[1,0]:(szy - offsets[1,1]):2 ]
        sig_out[:endx, :endy] = sig_in[ ix_(indsx, indsy) ]
    sig_in[:] = sig_out[:]
    
def _reorganize_coeffs3_forward(sig_in, m, offsets, data_layout):
    if data_layout.lower() == 'resolution':
        sig_out = zeros_like(sig_in)
        indsx = arange(shape(sig_in)[0])
        indsy = arange(shape(sig_in)[1])
        indsz = arange(shape(sig_in)[2])
        endx = shape(sig_in)[0]; endy = shape(sig_in)[1]; endz = shape(sig_in)[2]
        for res in range(1, m + 1):
            szx = len(indsx); szy = len(indsy); szz = len(indsz)
            psiinds_x = concatenate( ( indsx[:offsets[0,0]], indsx[(offsets[0,0] + 1):(szx - offsets[0,1]):2], indsx[(szx - offsets[0,1]):] ) ) # psi-indices
            psiinds_y = concatenate( ( indsy[:offsets[1,0]], indsy[(offsets[1,0] + 1):(szy - offsets[1,1]):2], indsy[(szy - offsets[1,1]):] ) )
            psiinds_z = concatenate( ( indsz[:offsets[2,0]], indsz[(offsets[2,0] + 1):(szz - offsets[2,1]):2], indsz[(szz - offsets[2,1]):] ) )
            phiinds_x = indsx[offsets[0,0]:(szx - offsets[0,1]):2]
            phiinds_y = indsy[offsets[1,0]:(szy - offsets[1,1]):2]
            phiinds_z = indsz[offsets[2,0]:(szz - offsets[2,1]):2]
            
            sig_out[ :(endx-len(psiinds_x)), :(endy-len(psiinds_y)), (endz-len(psiinds_z)):endz] =     sig_in[ ix_(phiinds_x, phiinds_y, psiinds_z) ]
            sig_out[ :(endx-len(psiinds_x)), (endy-len(psiinds_y)):endy, :(endz-len(psiinds_z))] =     sig_in[ ix_(phiinds_x, psiinds_y, phiinds_z) ]
            sig_out[ (endx-len(psiinds_x)):endx, :(endy-len(psiinds_y)), :(endz-len(psiinds_z))] =     sig_in[ ix_(psiinds_x, phiinds_y, phiinds_z) ]
            
            sig_out[ (endx-len(psiinds_x)):endx, (endy-len(psiinds_y)):endy, :(endz-len(psiinds_z))] =     sig_in[ ix_(psiinds_x, psiinds_y, phiinds_z) ]
            sig_out[ (endx-len(psiinds_x)):endx, :(endy-len(psiinds_y)), (endz-len(psiinds_z)):endz] =     sig_in[ ix_(psiinds_x, phiinds_y, psiinds_z) ]
            sig_out[ :(endx-len(psiinds_x)), (endy-len(psiinds_y)):endy, (endz-len(psiinds_z)):endz] =     sig_in[ ix_(phiinds_x, psiinds_y, psiinds_z) ]
            
            sig_out[ (endx-len(psiinds_x)):endx, (endy-len(psiinds_y)):endy, (endz-len(psiinds_z)):endz] = sig_in[ ix_(psiinds_x, psiinds_y, psiinds_z) ]
            
            endx -= len(psiinds_x); endy -= len(psiinds_y); endz -= len(psiinds_z)
            indsx = indsx[ offsets[0,0]:(szx - offsets[0,1]):2 ] 
            indsy = indsy[ offsets[1,0]:(szy - offsets[1,1]):2 ]
            indsz = indsz[ offsets[2,0]:(szz - offsets[2,1]):2 ]
        sig_out[:endx, :endy, :endz] = sig_in[ ix_(indsx, indsy, indsz) ]
    sig_in[:] = sig_out[:]

def _reorganize_coeffs_reverse(x, m, offsets, data_layout):
    inds = arange(shape(x)[0])
    y = zeros_like(x)
    resstart = array(zeros((1,m + 1), int))
    resend   = array(zeros((1,m + 1), int))
    resstart[0,0] = inds[0]; resend[0,0] = inds[-1]
    if data_layout.lower() == 'time':
        for res in range(1, m + 1):
            sz = len(inds)
            inds = inds[ offsets[0,0]:(sz - offsets[0,1]):2 ]
            resstart[0, res] = inds[0]
            resend[0, res] = inds[-1]
    if data_layout.lower() == 'resolution':
        endy = shape(x)[0]
        for res in range(1, m + 1):
            sz = len(inds)
            xindices = concatenate( ( inds[:offsets[0,0]], inds[(offsets[0,0] + 1):(sz - offsets[0,1]):2], inds[(sz - offsets[0,1]):] ) ) # psi-indices
            resstart[0,res] = inds[offsets[0,0]]; resend[0,res]   = inds[sz - offsets[0,1] - 1]
            y[xindices] = x[(endy-len(xindices)):endy]
            endy = endy-len(xindices)
            inds = inds[offsets[0,0]:(sz - offsets[0,1]):2]
        y[inds]= x[:endy]
        x[:] = y[:]
    return resstart, resend




     
def _reorganize_coeffs2_reverse(sig_in, m, offsets, data_layout):
    indsx = arange(shape(sig_in)[0])
    indsy = arange(shape(sig_in)[1])
    sig_out = zeros_like(sig_in)
    resstart = array(zeros((2,m + 1), int))
    resend   = array(zeros((2,m + 1), int))
    resstart[0,0] = indsx[0]; resend[0,0] = indsx[-1]
    resstart[1,0] = indsy[0]; resend[1,0] = indsy[-1]
    if data_layout.lower() == 'time':
        for res in range(1, m + 1):
            szx = len(indsx); szy = len(indsy)
            indsx = indsx[ offsets[0,0]:(szx - offsets[0,1]):2 ]
            indsy = indsy[ offsets[1,0]:(szy - offsets[1,1]):2 ]
            resstart[0,res] = indsx[0]; resend[0,res] = indsx[-1]
            resstart[1,res] = indsy[0]; resend[1,res] = indsy[-1]
    if data_layout.lower() == 'resolution':
        endx = shape(sig_in)[0]; endy = shape(sig_in)[1]
        for res in range(1, m + 1):
            szx = len(indsx); szy = len(indsy)
            psiinds_x = concatenate( ( indsx[:offsets[0,0]], indsx[(offsets[0,0] + 1):(szx - offsets[0,1]):2], indsx[(szx - offsets[0,1]):] ) ) # psi-indices
            psiinds_y = concatenate( ( indsy[:offsets[1,0]], indsy[(offsets[1,0] + 1):(szy - offsets[1,1]):2], indsy[(szy - offsets[1,1]):] ) )
            phiinds_x = indsx[offsets[0,0]:(szx - offsets[0,1]):2]
            phiinds_y = indsy[offsets[1,0]:(szy - offsets[1,1]):2]
            
            resstart[0, res] = indsx[offsets[0,0]]; resend[0, res]   = indsx[szx - offsets[0,1] - 1]
            resstart[1, res] = indsy[offsets[1,0]]; resend[1, res]   = indsy[szy - offsets[1,1] - 1]
            
            sig_out[ ix_(psiinds_x, phiinds_y) ] = sig_in[ (endx-len(psiinds_x)):endx, :(endy-len(psiinds_y))]     
            sig_out[ ix_(phiinds_x, psiinds_y) ] = sig_in[ :(endx-len(psiinds_x)), (endy-len(psiinds_y)):endy]
            sig_out[ ix_(psiinds_x, psiinds_y) ] = sig_in[ (endx-len(psiinds_x)):endx, (endy-len(psiinds_y)):endy]
            
            endx = endx - len(psiinds_x); endy = endy - len(psiinds_y)
            indsx = indsx[offsets[0,0]:(szx - offsets[0,1]):2]
            indsy = indsy[offsets[1,0]:(szy - offsets[1,1]):2]
        sig_out[ix_(indsx, indsy)] = sig_in[:endx, :endy]
        sig_in[:] = sig_out[:]
    return resstart, resend
            
def _reorganize_coeffs3_reverse(sig_in, m, offsets, data_layout):
    indsx = arange(shape(sig_in)[0])
    indsy = arange(shape(sig_in)[1])
    indsz = arange(shape(sig_in)[2])
    sig_out = zeros_like(sig_in)
    resstart = array(zeros((3,m + 1), int))
    resend   = array(zeros((3,m + 1), int))
    resstart[0,0] = indsx[0]; resend[0,0] = indsx[-1]
    resstart[1,0] = indsy[0]; resend[1,0] = indsy[-1]
    resstart[2,0] = indsz[0]; resend[2,0] = indsz[-1]
    if data_layout.lower() == 'time':
        for res in range(1, m + 1):
            szx = len(indsx); szy = len(indsy); szz = len(indsz)
            indsx = indsx[ offsets[0,0]:(szx - offsets[0,1]):2 ]
            indsy = indsy[ offsets[1,0]:(szy - offsets[1,1]):2 ]
            indsz = indsz[ offsets[2,0]:(szz - offsets[2,1]):2 ]
            resstart[0,res] = indsx[0]; resend[0,res] = indsx[-1]
            resstart[1,res] = indsy[0]; resend[1,res] = indsy[-1]
            resstart[2,res] = indsz[0]; resend[2,res] = indsz[-1]
    if data_layout.lower() == 'resolution':
        endx = shape(sig_in)[0]; endy = shape(sig_in)[1]; endz = shape(sig_in)[2]
        for res in range(1, m + 1):
            szx = len(indsx); szy = len(indsy); szz = len(indsz)
            psiinds_x = concatenate( ( indsx[:offsets[0,0]], indsx[(offsets[0,0] + 1):(szx - offsets[0,1]):2], indsx[(szx - offsets[0,1]):] ) ) # psi-indices
            psiinds_y = concatenate( ( indsy[:offsets[1,0]], indsy[(offsets[1,0] + 1):(szy - offsets[1,1]):2], indsy[(szy - offsets[1,1]):] ) )
            psiinds_z = concatenate( ( indsz[:offsets[2,0]], indsz[(offsets[2,0] + 1):(szz - offsets[2,1]):2], indsz[(szz - offsets[2,1]):] ) )
            phiinds_x = indsx[offsets[0,0]:(szx - offsets[0,1]):2]
            phiinds_y = indsy[offsets[1,0]:(szy - offsets[1,1]):2]
            phiinds_z = indsz[offsets[2,0]:(szz - offsets[2,1]):2]
            
            resstart[0, res] = indsx[offsets[0,0]]; resend[0, res]   = indsx[szx - offsets[0,1] - 1]
            resstart[1, res] = indsy[offsets[1,0]]; resend[1, res]   = indsy[szy - offsets[1,1] - 1]
            resstart[2, res] = indsz[offsets[2,0]]; resend[2, res]   = indsz[szz - offsets[2,1] - 1]
            
            sig_out[ ix_(phiinds_x, phiinds_y, psiinds_z) ] = sig_in[ :(endx-len(psiinds_x)), :(endy-len(psiinds_y)), (endz-len(psiinds_z)):endz]
            sig_out[ ix_(phiinds_x, psiinds_y, phiinds_z) ] = sig_in[ :(endx-len(psiinds_x)), (endy-len(psiinds_y)):endy, :(endz-len(psiinds_z))]
            sig_out[ ix_(psiinds_x, phiinds_y, phiinds_z) ] = sig_in[ (endx-len(psiinds_x)):endx, :(endy-len(psiinds_y)), :(endz-len(psiinds_z))]
            
            sig_out[ ix_(psiinds_x, psiinds_y, phiinds_z) ] = sig_in[ (endx-len(psiinds_x)):endx, (endy-len(psiinds_y)):endy, :(endz-len(psiinds_z))]
            sig_out[ ix_(psiinds_x, phiinds_y, psiinds_z) ] = sig_in[ (endx-len(psiinds_x)):endx, :(endy-len(psiinds_y)), (endz-len(psiinds_z)):endz]
            sig_out[ ix_(phiinds_x, psiinds_y, psiinds_z) ] = sig_in[ :(endx-len(psiinds_x)), (endy-len(psiinds_y)):endy, (endz-len(psiinds_z)):endz]
            
            sig_out[ ix_(psiinds_x, psiinds_y, psiinds_z) ] = sig_in[ (endx-len(psiinds_x)):endx, (endy-len(psiinds_y)):endy, (endz-len(psiinds_z)):endz]
            
            endx = endx - len(psiinds_x); endy = endy - len(psiinds_y); endz = endz - len(psiinds_z)
            indsx = indsx[offsets[0,0]:(szx - offsets[0,1]):2]
            indsy = indsy[offsets[1,0]:(szy - offsets[1,1]):2]
            indsz = indsz[offsets[2,0]:(szz - offsets[2,1]):2]
        sig_out[ix_(indsx, indsy, indsz)] = sig_in[:endx, :endy, :endz]
        sig_in[:] = sig_out[:]
    return resstart, resend
            
            
            
def tensor2_impl(x, fx, fy, bd_mode):
    sz = arange(len(shape(x)))
    fx(x, bd_mode)
    y = transpose(x, concatenate( ([1, 0], sz[2:]) ))
    fy(y, bd_mode)
    x[:] = transpose(y, concatenate( ([1, 0], sz[2:]) ))
# End tensor2_impl

def tensor3_impl(x, fx, fy, fz, bd_mode):
    sz = arange(len(shape(x)))
    fx(x, bd_mode)
    y =    transpose(x, concatenate( ([1, 2, 0], sz[3:]) ))
    fy(y, bd_mode)
    z =    transpose(y, concatenate( ([1, 2, 0], sz[3:]) ))
    fz(z, bd_mode)
    x[:] = transpose(z, concatenate( ([1, 2, 0], sz[3:]) ))
# End tensor3_impl

# testcode
  
    
def cascade_alg(m, a, b, wave_name, scaling, dual):
    coords = zeros((b-a)*2**m)
    if scaling:
        coords[0] = 1
    else:
        coords[b - a] = 1
    
    t = linspace(a, b, (b-a)*2**m)
    idwt_impl(coords, wave_name, m=m, bd_mode='per', dual=dual)
    coords = concatenate([coords[(b*2**m):((b-a)*2**m)], \
                          coords[0:(b*2**m)]])
    plt.plot(t, 2**(m/2.)*coords, 'k-')
    plt.show()
    plt.close()
# End cascade_alg

def freqresp_alg(wave_name, lowpass, dual):
    N = 128
    n = arange(N)
    omega = 2*pi*n/float(N)

    g = zeros(N)
    if lowpass:
        g[0] = 1
    else:
        g[1] = 1
    
    idwt_impl(g, wave_name, bd_mode='per', dual=dual, data_layout='time')
    plt.plot(omega, abs(fft.fft(g)), 'k-')
    plt.show()
    plt.close()    
# End freqresp_alg
    
def _test_dwt_different_sizes():
    for wave_name in ['cdf97', 'spline53', 'pwl0', 'pwl2', 'haar', 'spline4.4']:
        print('Testing the DWT on different input sizes')
        m = 4
    
        print('Testing 2D with one channel: %s' % wave_name)
        img = random.random((64,64))
        img2 = img.copy()
        dwt_impl(img2, wave_name, m=m, dims=2)
        idwt_impl(img2, wave_name, m=m, dims=2)
        # print(img2)
        diff = abs(img2-img).max()
        assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff
        
        print('Testing 2D with three channels: %s' % wave_name)
        img = random.random((64, 64, 3))
        img2 = img.copy()
        dwt_impl(img2, wave_name, m=m)
        idwt_impl(img2, wave_name, m=m)
        diff = abs(img2-img).max()
        assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff
        
        print('Testing 1D with one channel: %s' % wave_name)
        sd = random.random(64)
        sd2 = sd.copy()
        dwt_impl(sd2, wave_name, m=m)
        idwt_impl(sd2, wave_name, m=m)
        diff = abs(sd2-sd).max()
        assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff
        
        print('Testing 1D with two channels: %s' % wave_name)
        sd = random.random((64,2))
        sd2 = sd.copy()
        dwt_impl(sd2, wave_name, m=m)
        idwt_impl(sd2, wave_name, m=m)
        diff = abs(sd2-sd).max()
        assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff
        
        print('Testing 3D with one channel: %s' % wave_name)
        sd = random.random((64,64,64))
        sd2 = sd.copy()
        dwt_impl(sd2, wave_name, m=m, dims=3)
        idwt_impl(sd2, wave_name, m=m, dims=3)
        diff = abs(sd2-sd).max()
        assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff
        
        print('Testing 3D with two channels: %s' % wave_name)
        sd = random.random((64,64,64,3))
        sd2 = sd.copy()
        dwt_impl(sd2, wave_name, m=m)
        idwt_impl(sd2, wave_name, m=m)
        diff = abs(sd2-sd).max()
        assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff
    
    
    
def _test_orthogonality():
    for wave_name in ['db2', 'db4']:
        print('Testing orthonormal wavelets:')
        x0 = random.random(32)
        
        print('Testing that the IDWT inverts the DWT: %s' % wave_name)
        x = x0.copy()
        dwt_impl(x, wave_name, m=2, bd_mode='per')
        idwt_impl(x, wave_name, m=2, bd_mode='per')
        diff = abs(x-x0).max()
        assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff
    
        print('See that the wavelet transform equals the dual wavelet transform: %s' % wave_name)
        x = x0.copy()
        dwt_impl(x, wave_name, m=2, bd_mode='per', dual=True)
        dwt_impl(x0, wave_name, m=2, bd_mode='per', dual=False)
        diff = abs(x-x0).max()
        assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff
    
        print('Apply the transpose, to see that the transpose equals the inverse: %s' % wave_name)
        x = x0.copy()
        dwt_impl(x, wave_name, m=2, bd_mode='per', transpose=True)
        dwt_impl(x, wave_name, m=2, bd_mode='per', transpose=False)
        diff = abs(x-x0).max()
        assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff

def _test_simple_dwt2():
    print('Testing simple DWT2')
    img2 = random.random((32, 32, 3))
    img = img2.copy()
# Begin simple_dwt2
    f = lambda x, bd_mode: dwt_impl(x, 'cdf97', m=4, bd_mode=bd_mode, dims=1)
    tensor2_impl(img, f, f, 'symm')
# End simple_dwt2
# Begin simple_idwt2
    invf = lambda x, bd_mode: idwt_impl(x,'cdf97',m=4,bd_mode=bd_mode, dims=1)
    tensor2_impl(img, invf, invf, 'symm')
# End simple_idwt2
    diff = abs(img2-img).max()
    assert (diff < 1E-13 and diff != 0) , 'bug, diff=%s' % diff
    
if __name__=='__main__':
    _test_orthogonality()
    _test_dwt_different_sizes()
    _test_simple_dwt2()