import random
from numpy import *
import matplotlib.pyplot as plt
from sound import *
from images import *
import scipy.io as sio
import os
import os.path




def dwt_impl(x, nres, wave_name, forward=True, dim=1, bd_mode='symm', dual=False, transpose=False):
    # For classes
    
    f=find_kernel(wave_name, forward, dual, transpose)
    if transpose:
        forward = not forward
    if forward:
        if dim == 2:
            DWT2Impl_internal(x, nres, f, bd_mode)
        elif dim == 1:
            DWTImpl_internal(x, nres, f, bd_mode)
    else:
        if dim == 2:
            IDWT2Impl_internal(x, nres, f, bd_mode)
        elif dim == 1:
            IDWTImpl_internal(x, nres, f, bd_mode)

def DWT2Impl(x, nres, wave_name, bd_mode='symm', dual=False, transpose=False):
    """
    x:         Matrix whose DWT will be computed along the first dimension(s).      
    nres:      Number of resolutions.
    wave_name: Name of the wavelet. Possible names are:
               'cdf97' - CDF 9/7 wavelet
               'cdf53' - Spline 5/3 wavelet  
               'pwl0'  - Piecewise linear wavelets with 0 vanishing moments
               'pwl2'  - Piecewise linear wavelets with 2 vanishing moments
               'Haar'  - The Haar wavelet
               'dbX'   - Dauberchies orthnormal wavelet with X vanishing
                         moments
               'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
                         with X vanishing moments
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'bd'     - Boundary wavelets
               'bd_pre' - Boundary wavelets with preconditioning
    dual:      Whether to apply the dual wavelet rather than the wavelet itself. Default: False
    transpose: Whether the transpose is to be taken. Default: False
    """
    f = find_kernel(wave_name, True, dual, transpose)
    if transpose:
        IDWT2Impl_internal(x, nres, f, bd_mode)
    else:
        DWT2Impl_internal(x, nres, f, bd_mode)
        
def DWTImpl(x, nres, wave_name, bd_mode='symm', dual=False, transpose=False):
    """
    x:         Matrix whose DWT will be computed along the first dimension(s).      
    nres:      Number of resolutions.
    wave_name: Name of the wavelet. Possible names are:
               'cdf97' - CDF 9/7 wavelet
               'cdf53' - Spline 5/3 wavelet  
               'pwl0'  - Piecewise linear wavelets with 0 vanishing moments
               'pwl2'  - Piecewise linear wavelets with 2 vanishing moments
               'Haar'  - The Haar wavelet
               'dbX'   - Dauberchies orthnormal wavelet with X vanishing
                         moments
               'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
                         with X vanishing moments
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'bd'     - Boundary wavelets
               'bd_pre' - Boundary wavelets with preconditioning
    dual:      Whether to apply the dual wavelet rather than the wavelet itself. Default: False
    transpose: Whether the transpose is to be taken. Default: False
    """
    f = find_kernel(wave_name, True, dual, transpose)
    if transpose:
        IDWTImpl_internal(x, nres, f, bd_mode)
    else:
        DWTImpl_internal(x, nres, f, bd_mode)

def IDWT2Impl(x, nres, wave_name, bd_mode='symm', dual=False, transpose=False):
    """
    x:         Matrix whose IDWT will be computed along the first dimension(s).      
    nres:      Number of resolutions.
    wave_name: Name of the wavelet. Possible names are:
               'cdf97' - CDF 9/7 wavelet
               'cdf53' - Spline 5/3 wavelet  
               'pwl0'  - Piecewise linear wavelets with 0 vanishing moments
               'pwl2'  - Piecewise linear wavelets with 2 vanishing moments
               'Haar'  - The Haar wavelet
               'dbX'   - Dauberchies orthnormal wavelet with X vanishing
                         moments
               'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
                         with X vanishing moments
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'bd'     - Boundary wavelets
               'bd_pre' - Boundary wavelets with preconditioning
    dual:      Whether to apply the dual wavelet rather than the wavelet itself. Default: False
    transpose: Whether the transpose is to be taken. Default: False
    """
    f = find_kernel(wave_name, False, dual, transpose)
    if transpose:
        DWT2Impl_internal(x, nres, f, bd_mode)
    else:
        IDWT2Impl_internal(x, nres, f, bd_mode)
        
def IDWTImpl(x, nres, wave_name, bd_mode='symm', dual=False, transpose=False):
    """
    x:         Matrix whose IDWT will be computed along the first dimension(s).      
    nres:      Number of resolutions.
    wave_name: Name of the wavelet. Possible names are:
               'cdf97' - CDF 9/7 wavelet
               'cdf53' - Spline 5/3 wavelet  
               'pwl0'  - Piecewise linear wavelets with 0 vanishing moments
               'pwl2'  - Piecewise linear wavelets with 2 vanishing moments
               'Haar'  - The Haar wavelet
               'dbX'   - Dauberchies orthnormal wavelet with X vanishing
                         moments
               'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
                         with X vanishing moments
    bd_mode:   Boundary extension mode. Possible modes are. 
               'per'    - Periodic extension
               'symm'   - Symmetric extension (default)
               'bd'     - Boundary wavelets
               'bd_pre' - Boundary wavelets with preconditioning
    dual:      Whether to apply the dual wavelet rather than the wavelet itself. Default: False
    transpose: Whether the transpose is to be taken. Default: False
    """
    f = find_kernel(wave_name, False, dual, transpose)
    if transpose:
        DWTImpl_internal(x, nres, f, bd_mode)
    else:
        IDWTImpl_internal(x, nres, f, bd_mode)
        
        
        
            
def DWT2Impl_internal(X, nres, f, bd_mode):
    M, N = shape(X)[0:2]
    for res in range(nres):
        for n in range(0,N,2**res): 
            f(X[0::2**res, n], bd_mode)
        for m in range(0,M,2**res):
            f(X[m, 0::2**res], bd_mode)
    reorganize_coeffs2_forward(X, nres)   
            
def IDWT2Impl_internal(X, nres, f, bd_mode):
    reorganize_coeffs2_reverse(X, nres)   
    M, N = shape(X)[0:2]        
    for res in range(nres - 1, -1, -1):
        for n in range(0, N, 2**res):
            f(X[0::2**res, n], bd_mode)
        for m in range(0, M, 2**res):
            f(X[m, 0::2**res], bd_mode)
  
def DWTImpl_internal(x, m, f, bd_mode):
    for res in range(m):
        f(x[0::2**res], bd_mode)
    reorganize_coeffs_forward(x, m)
            
def IDWTImpl_internal(x, m, f, bd_mode):
    reorganize_coeffs_reverse(x, m)
    for res in range(m - 1, -1, -1):
        f(x[0::2**res], bd_mode)
        
        
        
        
        
def find_kernel(wave_name, forward, dual, transpose):
    if transpose:
        forward = not forward
        dual = not dual
    if forward:
        if dual:
            f = find_kernel_dwt_dual(wave_name)
        else:
            f = find_kernel_dwt(wave_name)
    else:
        if dual:
            f = find_kernel_idwt_dual(wave_name)
        else:
            f = find_kernel_idwt(wave_name)
    return f



def find_kernel_dwt_dual(wave_name):
    f = 0
    if wave_name.lower() =='cdf97':
        f = dwt_kernel_97_dual
    elif wave_name.lower() == 'cdf53':
        f = dwt_kernel_53_dual
    elif wave_name.lower() == 'pwl0':
        f = dwt_kernel_pwl0_dual
    elif wave_name.lower() == 'pwl2':
        f = dwt_kernel_pwl2_dual
    elif wave_name.lower() == 'haar':
        f = dwt_kernel_haar
    elif wave_name[:2].lower() == 'db' and not wave_name[-1].lower() =='x':
        vm = int(wave_name[2::])
        filters = getDBfilter(vm, 0)
        f = lambda x, bd_mode: dwt_kernel_ortho_dual(x, filters, bd_mode)
    elif wave_name[:2].lower() == 'db':
        vm = int(wave_name[2:-1])
        filters = liftingfactortho(vm, 0, 1)
        f = lambda x, bd_mode: dwt_kernel_ortho_dual(x, filters, bd_mode)
    elif wave_name[:3].lower() == 'sym' and not wave_name[-1].lower() =='x':
        vm = int(wave_name[3::])
        filters = getDBfilter(vm, 1)
        f = lambda x, bd_mode: det_kernel_ortho_dual(x, filters, bd_mode)
    elif wave_name[:3].lower() == 'sym':
        vm = int(wave_name[3:-1])
        filters = liftingfactortho(vm, 1, 1)
        f = lambda x, bd_mode: dwt_kernel_ortho_dual(x, filters, bd_mode)
    return f

def find_kernel_dwt(wave_name):
    f = 0
    if wave_name.lower() =='cdf97':
        f = dwt_kernel_97
    elif wave_name.lower() == 'cdf53':
        f = dwt_kernel_53
    elif wave_name.lower() == 'pwl0':
        f = dwt_kernel_pwl0
    elif wave_name.lower() == 'pwl2':
        f = dwt_kernel_pwl2
    elif wave_name.lower() == 'haar':
        f = dwt_kernel_haar
    elif wave_name[:2].lower() == 'db' and wave_name[-1].lower() != 'x':
        vm = int(wave_name[2::])
        filters = getDBfilter(vm, 0)
        f = lambda x, bd_mode: dwt_kernel_ortho(x, filters, bd_mode)
    elif wave_name[:2].lower() == 'db':
        vm = int(wave_name[2:-1])
        filters = liftingfactortho(vm, 0, 1)
        f = lambda x, bd_mode: dwt_kernel_ortho(x, filters, bd_mode)
    elif wave_name[:3].lower() == 'sym' and not wave_name[-1].lower() =='x':
        vm = int(wave_name[3:])
        filters = getDBfilter(vm, 1)
        f = lambda x, bd_mode: dwt_kernel_ortho(x, filters, bd_mode)
    elif wave_name[:3].lower() == 'sym':
        vm = int(wave_name[3:-1])
        filters = liftingfactortho(vm, 1, 1)
        f = lambda x, bd_mode: dwt_kernel_ortho(x, filters, bd_mode)
    return f
    
def find_kernel_idwt_dual(wave_name):
    f = 0
    if wave_name.lower() =='cdf97':
        f = idwt_kernel_97_dual
    elif wave_name.lower() == 'cdf53':
        f = idwt_kernel_53_dual
    elif wave_name.lower() == 'pwl0':
        f = idwt_kernel_pwl0_dual
    elif wave_name.lower() == 'pwl2':
        f = idwt_kernel_pwl2_dual
    elif wave_name.lower() == 'haar':
        f = idwt_kernel_haar
    elif wave_name[:2].lower() == 'db' and not wave_name[-1].lower() =='x':
        vm = int(wave_name[2::])
        filters = getDBfilter(vm, 0)
        f = lambda x, bd_mode: idwt_kernel_ortho_dual(x, filters, bd_mode)
    elif wave_name[:2].lower() == 'db':
        vm = int(wave_name[2:-1])
        filters = liftingfactortho(vm, 0, 1)
        f = lambda x, bd_mode: idwt_kernel_ortho_dual(x, filters, bd_mode)
    elif wave_name[:3].lower() == 'sym' and not wave_name[-1].lower() =='x':
        vm = int(wave_name[3::])
        filters = getDBfilter(vm, 1)
        f = lambda x, bd_mode: idwt_kernel_ortho_dual(x, filters, bd_mode)
    elif wave_name[:3].lower() == 'sym':
        vm = int(wave_name[3:-1])
        filters = liftingfactortho(vm, 1, 1)
        f = lambda x, bd_mode: idwt_kernel_ortho_dual(x, filters, bd_mode)
    return f

def find_kernel_idwt(wave_name):
    f = 0
    if wave_name.lower() =='cdf97':
        f = idwt_kernel_97
    elif wave_name.lower() == 'cdf53':
        f = idwt_kernel_53
    elif wave_name.lower() == 'pwl0':
        f = idwt_kernel_pwl0
    elif wave_name.lower() == 'pwl2':
        f = idwt_kernel_pwl2
    elif wave_name.lower() == 'haar':
        f = idwt_kernel_haar
    elif wave_name[:2].lower() == 'db' and not wave_name[-1].lower() =='x':
        vm = int(wave_name[2::])
        filters = getDBfilter(vm, 0)
        f = lambda x, bd_mode: idwt_kernel_ortho(x, filters, bd_mode)
    elif wave_name[:2].lower() == 'db':
        vm = int(wave_name[2:-1])
        filters = liftingfactortho(vm, 0, 1)
        f = lambda x, bd_mode: idwt_kernel_ortho(x, filters, bd_mode)
    elif wave_name[:3].lower() == 'sym' and not wave_name[-1].lower() =='x':
        vm = int(wave_name[3::])
        filters = getDBfilter(vm, 1)
        f = lambda x, bd_mode: idwt_kernel_ortho(x, filters, bd_mode)
    elif wave_name[:3].lower() == 'sym':
        vm = int(wave_name[3:-1])
        filters = liftingfactortho(vm, 1, 1)
        f = lambda x, bd_mode: idwt_kernel_ortho(x, filters, bd_mode)
    return f

def getDBfilter(vm, type):
    filter = 0;
    dest = '../matlab/var' 
    if (type == 0):
        filename = '%s/DB%d.mat' % (dest, vm)
    else:
        filename = '%s/sym%d.mat' % (dest, vm)
    if os.path.isfile(filename):
        vals = sio.loadmat(filename)
        newfilter = vals['filter']
        filter = {'lambdas': newfilter['lambdas'][0][0],'alpha': float(newfilter['alpha']),'beta': float(newfilter['beta'])}
    else:
        filter = liftingfactortho(vm, type)
        
        if not os.path.isdir(dest):
            os.mkdir(dest)
        
        sio.savemat(filename, {'filter': filter})
    return filter

    
def dwt_kernel_filters(H0, H1, G0, G1, x, bd_mode):
    symm = bd_mode.lower() =='symm'
    f0, f1 = H0, H1
    x0 = x.copy()
    x1 = x.copy()
    filterS(f0, x0, symm)
    filterS(f1, x1, symm)
    x[::2] = x0[::2]
    x[1::2] = x1[1::2]
    
def dwt_kernel_filters_dual(H0, H1, G0, G1, x, bd_mode):
    symm = bd_mode.lower() =='symm'
    f0, f1 = G0, G1
    x0 = x.copy()
    x1 = x.copy()
    filterS(f0, x0, symm)
    filterS(f1, x1, symm)
    x[::2] = x0[::2]
    x[1::2] = x1[1::2]

def idwt_kernel_filters(H0, H1, G0, G1, x, bd_mode):
    symm = bd_mode.lower() =='symm'
    f0, f1 = G0, G1
    x0 = x.copy(); x0[1::2] = 0
    x1 = x.copy(); x1[::2] = 0
    filterS(f0, x0, symm)
    filterS(f1, x1, symm)
    x[:] = x0 + x1
    
def idwt_kernel_filters_dual(H0, H1, G0, G1, x, bd_mode):
    symm = bd_mode.lower() =='symm'
    f0, f1 = H0, H1
    x0 = x.copy(); x0[1::2] = 0
    x1 = x.copy(); x1[::2] = 0
    filterS(f0, x0, symm)
    filterS(f1, x1, symm)
    x[:] = x0 + x1
    
    
        
def reorganize_coeffs_forward(x, nres):
    N = shape(x)[0]
    y = zeros_like(x)
    sz = shape(x[0::2**nres])[0]
    y[0:sz] = x[0::2**nres]
    for res in range(nres, 0, -1):
        lw = shape(x[2**(res - 1)::2**res])[0]
        y[sz:(sz + lw)] = x[2**(res - 1)::2**res]
        sz += lw
    x[:] = y[:]
    
def reorganize_coeffs_reverse(x, nres):
    N = shape(x)[0]
    y = zeros_like(x)
    sz = shape(x[0::2**nres])[0]
    y[0::2**nres] = x[0:sz]
    for res in range(nres, 0, -1):
        lw = shape(x[2**(res - 1)::2**res])[0]
        y[2**(res - 1)::2**res] = x[sz:(sz + lw)]
        sz += lw
    x[:] = y[:]
    
def reorganize_coeffs2_forward(X, nres):
    M, N = shape(X)[0:2]
    Y = zeros_like(X)
    lc1, lc2 = shape(X[0::2**nres, 0::2**nres])[0:2]
    Y[0:lc1, 0:lc2] = X[0::2**nres, 0::2**nres]
    for res in range(nres, 0, -1):
        lw1, lw2 = shape(X[2**(res - 1)::2**res, 2**(res - 1)::2**res])[0:2]
        Y[lc1:(lc1 + lw1), 0:lc2] = X[2**(res - 1)::2**res, 0::2**res]
        Y[lc1:(lc1 + lw1), lc2:(lc2 + lw2)] = X[2**(res - 1)::2**res, 2**(res - 1)::2**res]
        Y[0:lc1, lc2:(lc2 + lw2)] = X[0::2**res, 2**(res - 1)::2**res]
        lc1 += lw1
        lc2 += lw2
    X[:] = Y[:]
    
def reorganize_coeffs2_reverse(X, nres):
    M, N = shape(X)[0:2]
    Y = zeros_like(X)
    lc1, lc2 = shape(X[0::2**nres, 0::2**nres])[0:2]
    Y[0::2**nres, 0::2**nres] = X[0:lc1, 0:lc2]
    for res in range(nres, 0, -1):
        lw1, lw2 = shape(X[2**(res - 1)::2**res, 2**(res - 1)::2**res])[0:2]
        Y[2**(res - 1)::2**res, 0::2**res] = X[lc1:(lc1 + lw1), 0:lc2]
        Y[2**(res - 1)::2**res, 2**(res - 1)::2**res] = X[lc1:(lc1 + lw1), lc2:(lc2 + lw2)]
        Y[0::2**res, 2**(res - 1)::2**res] = X[0:lc1, lc2:(lc2 + lw2)]
        lc1 += lw1
        lc2 += lw2
    X[:] = Y[:]
      
# Generic DWT/IDWT implementations


            
# Lifting steps
            
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
  
def lifting_odd_symm(lmbda, x, bd_mode):
    if (bd_mode.lower() == 'per') and mod(len(x), 2)!=0:
        raise AssertionError()
    x[1:-1:2] += lmbda*(x[0:-2:2] + x[2::2])
    if mod(len(x), 2)==0:
        if bd_mode.lower() == 'symm':
            x[-1] += 2*lmbda*x[-2] # With symmetric extension
        else:
            x[-1] += lmbda*(x[0]+x[-2])

def lifting_even(lmbda1, lmbda2, x):
    if mod(len(x), 2)!=0:
        raise AssertionError()
    x[0] += lmbda1*x[1] + lmbda2*x[-1]
    x[2:-1:2] += lmbda1*x[3::2] + lmbda2*x[1:-2:2]
            
def lifting_odd(lmbda1, lmbda2, x):
    if mod(len(x), 2)!=0:
        raise AssertionError()
    x[1:-2:2] += lmbda1*x[2:-1:2] + lmbda2*x[0:-3:2]
    x[-1] += lmbda1*x[0] + lmbda2*x[-2]                                                

# The Haar wavelet

def dwt_kernel_haar(x, bd_mode):
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
         
def idwt_kernel_haar(x, bd_mode):
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

def dwt_kernel_pwl0_dual(x, bd_mode):
    x /= sqrt(2)
    lifting_even_symm(0.5, x, bd_mode)
        
def dwt_kernel_pwl0(x, bd_mode):
    x *= sqrt(2)
    lifting_odd_symm(-0.5, x, bd_mode)
        
def idwt_kernel_pwl0_dual(x, bd_mode):
    x *= sqrt(2)
    lifting_even_symm(-0.5, x, bd_mode)
        
def idwt_kernel_pwl0(x, bd_mode):
    x /= sqrt(2)
    lifting_odd_symm(0.5, x, bd_mode)

def dwt_kernel_pwl2_dual(x, bd_mode):
    lifting_even_symm(0.5, x, bd_mode)
    lifting_odd_symm(-0.25, x, bd_mode)
    x /= sqrt(2)
    
def dwt_kernel_pwl2(x, bd_mode):
    lifting_odd_symm(-0.5, x, bd_mode)
    lifting_even_symm(0.25, x, bd_mode)
    x *= sqrt(2)
    
def idwt_kernel_pwl2_dual(x, bd_mode):
    x *= sqrt(2)
    lifting_odd_symm(0.25, x, bd_mode)
    lifting_even_symm(-0.5, x, bd_mode)
    
def idwt_kernel_pwl2(x, bd_mode):
    x /= sqrt(2)
    lifting_even_symm(-0.25, x, bd_mode)
    lifting_odd_symm(0.5, x, bd_mode)       
                
        
# JPEG2000-related wavelet kernels
        
def dwt_kernel_53_dual(x, bd_mode):
    x[0::2] *= 0.5
    x[1::2] *= 2
    lifting_even_symm(0.125, x, bd_mode)
    lifting_odd_symm(-1, x, bd_mode)
    
def dwt_kernel_53(x, bd_mode):
    x[0::2] *= 2
    x[1::2] *= 0.5
    lifting_odd_symm(-0.125, x, bd_mode)
    lifting_even_symm(1, x, bd_mode)
            
def idwt_kernel_53_dual(x, bd_mode):
    lifting_odd_symm(1, x, bd_mode)
    lifting_even_symm(-0.125, x, bd_mode)     
    x[0::2] *= 2
    x[1::2] *= 0.5
    
def idwt_kernel_53(x, bd_mode):
    lifting_even_symm(-1, x, bd_mode)
    lifting_odd_symm(0.125, x, bd_mode)     
    x[0::2] *= 0.5
    x[1::2] *= 2


def dwt_kernel_97_dual(x, bd_mode):
    lambda1=-0.586134342059950
    lambda2=-0.668067171029734
    lambda3=0.070018009414994
    lambda4=1.200171016244178
    alpha=-1.149604398860250
    beta=-0.869864451624777
    
    x[0::2] /= alpha
    x[1::2] /= beta
    lifting_even_symm(lambda4, x, bd_mode)
    lifting_odd_symm(lambda3, x, bd_mode)
    lifting_even_symm(lambda2, x, bd_mode)
    lifting_odd_symm(lambda1, x, bd_mode)
        
def dwt_kernel_97(x, bd_mode):
    lambda1=-0.586134342059950
    lambda2=-0.668067171029734
    lambda3=0.070018009414994
    lambda4=1.200171016244178
    alpha=-1.149604398860250
    beta=-0.869864451624777
    
    x[0::2] *= alpha
    x[1::2] *= beta
    lifting_odd_symm(-lambda4, x, bd_mode)
    lifting_even_symm(-lambda3, x, bd_mode)
    lifting_odd_symm(-lambda2, x, bd_mode)
    lifting_even_symm(-lambda1, x, bd_mode)
                
def idwt_kernel_97_dual(x, bd_mode):
    lambda1=-0.586134342059950
    lambda2=-0.668067171029734
    lambda3=0.070018009414994
    lambda4=1.200171016244178
    alpha=-1.149604398860250
    beta=-0.869864451624777
    
    lifting_odd_symm(-lambda1, x, bd_mode)
    lifting_even_symm(-lambda2, x, bd_mode)   
    lifting_odd_symm(-lambda3, x, bd_mode)
    lifting_even_symm(-lambda4, x, bd_mode)      
    x[0::2] *= alpha
    x[1::2] *= beta

def idwt_kernel_97(x, bd_mode):
    lambda1=-0.586134342059950
    lambda2=-0.668067171029734
    lambda3=0.070018009414994
    lambda4=1.200171016244178
    alpha=-1.149604398860250
    beta=-0.869864451624777
        
    lifting_even_symm(lambda1, x, bd_mode)
    lifting_odd_symm(lambda2, x, bd_mode)   
    lifting_even_symm(lambda3, x, bd_mode)
    lifting_odd_symm(lambda4, x, bd_mode)      
    x[0::2] /= alpha
    x[1::2] /= beta
        
# Orthonormal wavelets
        
def dwt_kernel_ortho_dual( x, filters, bd_mode):
    x[0::2] /= filters['alpha']
    x[1::2] /= filters['beta']
    for stepnr in range(filters['lambdas'].shape[0] - 1, 0, -2):
        lifting_odd(filters['lambdas'][stepnr, 1], filters['lambdas'][stepnr, 0], x)
        lifting_even(filters['lambdas'][stepnr -1, 1], filters['lambdas'][stepnr - 1, 0], x)   
    if mod(filters['lambdas'].shape[0], 2)==1:
        lifting_odd(filters['lambdas'][0, 1], filters['lambdas'][0, 0], x)

def dwt_kernel_ortho( x, filters, bd_mode):
    x[0::2] *= filters['alpha']
    x[1::2] *= filters['beta']
    for stepnr in range(filters['lambdas'].shape[0] - 1, 0, -2):
        lifting_even(-filters['lambdas'][stepnr, 0], -filters['lambdas'][stepnr, 1], x)
        lifting_odd(-filters['lambdas'][stepnr - 1, 0], -filters['lambdas'][stepnr - 1, 1], x)  
    if mod(filters['lambdas'].shape[0], 2)==1:
        lifting_even(-filters['lambdas'][0, 0], -filters['lambdas'][0, 1], x)
  
def idwt_kernel_ortho_dual( x, filters, bd_mode):
    stepnr = 0
    if mod(filters['lambdas'].shape[0], 2) == 1: # Start with an odd step
        lifting_odd(-filters['lambdas'][stepnr, 1], -filters['lambdas'][stepnr, 0], x)
        stepnr += 1
    while stepnr < filters['lambdas'].shape[0]:
        lifting_even(-filters['lambdas'][stepnr, 1], -filters['lambdas'][stepnr, 0], x)
        lifting_odd(-filters['lambdas'][stepnr + 1, 1], -filters['lambdas'][stepnr + 1, 0], x)
        stepnr += 2
    x[0::2] *= filters['alpha']
    x[1::2] *= filters['beta']
    
def idwt_kernel_ortho( x, filters, bd_mode):
    stepnr = 0
    if mod(filters['lambdas'].shape[0],2) == 1: # Start with an even step
        lifting_even(filters['lambdas'][stepnr, 0], filters['lambdas'][stepnr, 1], x)
        stepnr += 1
    while stepnr < filters['lambdas'].shape[0]:
        lifting_odd(filters['lambdas'][stepnr, 0], filters['lambdas'][stepnr, 1], x)
        lifting_even(filters['lambdas'][stepnr + 1, 0], filters['lambdas'][stepnr + 1, 1], x)
        stepnr += 2
    x[0::2] /= filters['alpha']
    x[1::2] /=filters['beta']
            



# testcode




def h0h1computeortho(N):
    vals = computeQN(N)
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
    return h0, h1

def liftingfactortho(N, type=0, debug_mode=False):
    """
    Assume that len(h1)==len(h0), and that h0 and h1 are even length and as symmetric as possible, with h0 with a minimum possible overweight of filter coefficients to the left, h1 to the right
    This function computes lifting steps l1, l2,...,ln, and constants alpha, beta so that ln ... l2 l1 H =  diag(alpha,beta), and stores these in files.
    This gives the following recipes for 
        Computing H: first multiply with diag(alpha,beta), then the inverses of the lifting steps in reverse order 
        Computing G: apply the lifting steps in the order given, finally multiply with diag(1/alpha,1/beta)
    ln is always odd, so that l1 is odd if and only if n is odd.
    All even lifting steps have only filter coefficients 0,1. All odd lifting steps have only filter coefficients -1,0
    """
    h0, h1 = h0h1computeortho(N)
    stepnr=0
    start1, end1, len1, start2, end2, len2 = 0, len(h0)/2-1, len(h0)/2,  0, len(h1)/2-1, len(h1)/2
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
  
    #[h00 h01; h10 h11]
    #convolve(h00,h11)-convolve(h10,h01)
    stepnr=stepnr+1

    # print [h00 h01; h10 h11], convolve(h00,h11)-convolve(h10,h01)
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
    
    # print [h00 h01; h10 h11], convolve(h00,h11)-convolve(h10,h01)
  
    # Add the final lifting, and compute alpha,beta
    alpha=sum(h00)
    beta=sum(h11)
    lastlift=-sum(h01)/beta
    if mod(len(h0)/2,2)==0:
        lambdas[stepnr,:] = [0,lastlift]
    else:
        lambdas[stepnr,:] = [lastlift,0]
    # [h00 h01; h10 h11]
    return {'lambdas': lambdas, 'alpha': alpha, 'beta': beta}
    
def computeQN(N):
    """
    Compute the coefficients of the polynomial Q^(N)((1-cos(w))/2).
    """
    QN=zeros(N)
    for k in range(N):
        QN[k] = 2*math.factorial(N+k-1)/(math.factorial(k)*math.factorial(N-1))
    vals = array([QN[0]])
    start = array([1.0])
    for k in range(1,N):
        start = convolve(start,[-1/4.0,1/2.0,-1/4.0])
        vals = hstack([0,vals])
        vals = hstack([vals,0])
        vals = vals + QN[k]*start
    return vals
    
    
    
def h0h1compute97():
    QN = computeQN(4)
    
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
    N= g0.shape[0]
    h1=g0*(-1)**(array(range(-(N-1)/2,(N+1)/2)))
    N= h0.shape[0]
    g1=h0*(-1)**(array(range(-(N-1)/2,(N+1)/2)))
    #print h0, h1, g0, g1
    return h0, h1
  
def liftingfact97():
    h0, h1 = h0h1compute97() # Should have 9 and 7 filter coefficients.
    h00, h01 = h0[0:9:2], h0[1:9:2]
    h10, h11 = h1[0:7:2], h1[1:7:2]
    
    lambdas=zeros(4)
    
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
    return lambdas, alpha, beta    
    
def _test_kernel(wave_name):
    print 'Testing %s, 1D' % wave_name
    res = random.random(16)
    x = zeros(16)
    x[:] = res[:]
    DWTImpl(x,2,wave_name)
    IDWTImpl(x,2,wave_name)
    diff = abs(x-res).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing %s, 2D' % wave_name
    res = random.random((16,2))
    x = zeros((16,2))
    x[:] = res[:]
    DWTImpl(x,2,wave_name)
    IDWTImpl(x,2,wave_name)
    diff = abs(x-res).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
def _test_kernel_ortho():
    print 'Testing orthonormal wavelets'
    res = random.random(16) # only this assumes that N is even
    x = zeros(16)
    
    print 'Testing that the reverse inverts the forward transform'
    x[0:16] = res[0:16]
    DWTImpl(x, 2, 'db4')
    IDWTImpl(x, 2, 'db4')
    diff = max(abs(x-res))
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing that the transform is orthogonal, i.e. that the transform and its dual are equal'
    x[0:16] = res[0:16]
    DWTImpl(x, 2, 'db4', 'per')
    DWTImpl(res, 2, 'db4', 'per', True)
    diff = max(abs(x-res))
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
def _test_dwt_different_sizes():
    print 'Testing the DWT on different input sizes'
    m = 4

    print 'Testing the DWT for greyscale image'
    img = random.random((32,32))
    img2 = zeros_like(img)
    img2[:] = img[:]
    DWT2Impl(img2, m, 'cdf97')
    IDWT2Impl(img2, m, 'cdf97')
    diff = abs(img2-img).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing the DWT for RGB image'
    img = random.random((32, 32, 3))
    img2 = zeros_like(img)
    img2[:] = img[:]
    DWT2Impl(img2, m, 'cdf97')
    IDWT2Impl(img2, m, 'cdf97')
    diff = abs(img2-img).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing the DWT for sound with one channel'
    sd = random.random(32)
    sd2 = zeros_like(sd)
    sd2[:] = sd[:]
    DWTImpl(sd2, m, 'cdf97')
    IDWTImpl(sd2, m, 'cdf97')
    diff = abs(sd2-sd).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing the DWT for sound with two channels'
    sd = random.random((32,2))
    sd2 = zeros_like(sd)
    sd2[:] = sd[:]
    DWTImpl(sd2, m, 'cdf97')
    IDWTImpl(sd2, m, 'cdf97')
    diff = abs(sd2-sd).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
def _test_orthogonality():
    print 'Testing that the wavelet and the dual wavelet are equal for orthonormal wavelets'
    x0 = random.random(32)
    
    print 'Testing that the IDWT inverts the DWT'
    x = x0.copy()
    DWTImpl(x, 2, 'db4', 'per', False)
    IDWTImpl(x, 2, 'db4', 'per', False)
    diff = abs(x-x0).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Apply the transpose, to see that the transpose equals the inverse'
    x = x0.copy()
    DWTImpl(x, 2, 'db4', 'per', False)
    IDWTImpl(x, 2, 'db4', 'per', True)
    diff = abs(x-x0).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff

    print 'See that the wavelet transform equals the dual wavelet transform'
    x = x0.copy()
    DWTImpl(x, 2, 'db4', 'per', True)
    DWTImpl(x0, 2, 'db4', 'per', False)
    diff = abs(x-x0).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff



if __name__=='__main__':
    _test_dwt_different_sizes()
    _test_kernel_ortho()
    _test_orthogonality()
    _test_kernel( 'cdf97')
    _test_kernel( 'cdf53')
    _test_kernel( 'pwl0')
    _test_kernel( 'pwl2')
    _test_kernel( 'haar')