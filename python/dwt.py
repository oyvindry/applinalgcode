import random
from numpy import *
import matplotlib.pyplot as plt
from sound import *
from images import *

def DWTKernelFilters(H0, H1, G0, G1, x, symm, dual):
    f0, f1 = H0, H1
    if dual:
        f0, f1 = G0, G1
    N = len(x)
    x0 = x.copy()
    x1 = x.copy()
    filterS(f0, x0, symm)
    filterS(f1, x1, symm)
    x[::2] = x0[::2]
    x[1::2] = x1[1::2]

def IDWTKernelFilters(H0, H1, G0, G1, x, symm, dual):
    f0, f1 = G0, G1
    if dual:
        f0, f1 = H0, H1
    N = len(x)
    x0 = x.copy(); x0[1::2] = 0
    x1 = x.copy(); x1[::2] = 0
    filterS(f0, x0, symm)
    filterS(f1, x1, symm)
    x[:] = x0 + x1
        
def reorganize_coefficients(x, nres, forward):
    """
    Permute one-dimensional DWT coefficients from the order used in-place DWT, to the order with low resolution coordinates first, as required by the DWT. If forward is false, the opposite reorganization is used. 
    
    x: The vector which holds the DWT coefficients. The reorganization is formed in the first dimension, but x may have a second dimension, as is the case for sound with more than one channel. 
    nres: The number of resolutions in x.
    forward: If forward is True, reorganize from in-place order to DWT coefficient order. If forward is False, reorganize from DWT coefficient order to in-place order.
    """
    N = shape(x)[0]
    y = zeros_like(x)
    sz = shape(x[0::2**nres])[0]
    if forward:
        y[0:sz] = x[0::2**nres]
    else:
        y[0::2**nres] = x[0:sz]
    for res in range(nres, 0, -1):
        lw = shape(x[2**(res - 1)::2**res])[0]
        if forward:
            y[sz:(sz + lw)] = x[2**(res - 1)::2**res]
        else:
            y[2**(res - 1)::2**res] = x[sz:(sz + lw)]
        sz += lw
    x[:] = y[:]
    
def reorganize_coefficients2(X, nres, forward):
    """
    Permute two-dimensional DWT coefficients from the order used in-place DWT, to the order with low resolution coordinates first, as required by the DWT. If forward is false, the opposite reorganization is used. 
    
    X: The matrix which holds the DWT coefficients. The reorganization is formed in the two first dimensions, but X may have a third dimension, as is the case for images with more than one colour component. 
    nres: The number of resolutions in X.
    forward: If forward is True, reorganize from in-place order to DWT coefficient order. If forward is False, reorganize from DWT coefficient order to in-place order.
    """
    M, N = shape(X)[0:2]
    Y = zeros_like(X)
    lc1, lc2 = shape(X[0::2**nres, 0::2**nres])[0:2]
    if forward:
        Y[0:lc1, 0:lc2] = X[0::2**nres, 0::2**nres]
    else:
        Y[0::2**nres, 0::2**nres] = X[0:lc1, 0:lc2]
    for res in range(nres, 0, -1):
        lw1, lw2 = shape(X[2**(res - 1)::2**res, 2**(res - 1)::2**res])[0:2]
        if forward:
            Y[lc1:(lc1 + lw1), 0:lc2] = X[2**(res - 1)::2**res, 0::2**res]
            Y[lc1:(lc1 + lw1), lc2:(lc2 + lw2)] = X[2**(res - 1)::2**res, 2**(res - 1)::2**res]
            Y[0:lc1, lc2:(lc2 + lw2)] = X[0::2**res, 2**(res - 1)::2**res]
        else:
            Y[2**(res - 1)::2**res, 0::2**res] = X[lc1:(lc1 + lw1), 0:lc2]
            Y[2**(res - 1)::2**res, 2**(res - 1)::2**res] = X[lc1:(lc1 + lw1), lc2:(lc2 + lw2)]
            Y[0::2**res, 2**(res - 1)::2**res] = X[0:lc1, lc2:(lc2 + lw2)]
        lc1 += lw1
        lc2 += lw2
    X[:] = Y[:]
      
# Generic DWT/IDWT implementations

def DWT2Impl(X, nres, f, symm=True, dual=False):
    """
    Compute a 2-dimensional DWT. The one-dimensional DWT is applied 
    to each row and column in X at each stage. X may have a third 
    axis, as is the case for images with more than one color 
    component. The DWT2 is then applied to each color component.
    
    X: A 2-dimensional object for which we apply the 2-dim DWT
    nres: The number of stages
    f: Wavelet kernel to apply. See DWTImpl for documentation
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    M, N = shape(X)[0:2]
    for res in range(nres):
        for n in range(0,N,2**res): 
            f(X[0::2**res, n], symm, dual)
        for m in range(0,M,2**res):
            f(X[m, 0::2**res], symm, dual)
    reorganize_coefficients2(X, nres, True)   
            
def IDWT2Impl(X, nres, f, symm=True, dual=False):
    """
    Compute a 2-dimensional IDWT. The one-dimensional IDWT is applied 
    to each row and column in X at each stage. X may have a third 
    axis, as is the case for images with more than one color 
    component. The IDWT2 is then applied to each color component.
    
    X: A 2-dimensional object for which we apply the 2-dim IDWT
    nres: The number of stages
    f: Wavelet kernel to apply. See IDWTImpl for documentation
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    reorganize_coefficients2(X, nres, False)   
    M, N = shape(X)[0:2]        
    for res in range(nres - 1, -1, -1):
        for n in range(0, N, 2**res):
            f(X[0::2**res, n], symm, dual)
        for m in range(0, M, 2**res):
            f(X[m, 0::2**res], symm, dual)
  
def DWTImpl(x, nres, f, symm=True, dual=False):
    """
    Compute the DWT of x for a given number of resolutions, using 
    wavelet kernel f. The kernel is assumed to compute one level 
    of the DWT in-place. DWTImpl is responsible for reorganizing the 
    output so that the low resolution coefficients comes first, as 
    required by the DWT. The DWT is computed along the first axis.
    x may have a second axis, as is the case for sound with more 
    than one channel. The DWT is then applied to each channel.
    
    x: The vector which we apply the DWT to.
    nres: The number of stages
    f: The wavelet kernel to apply. Supported kernels are 
        DWTKernelHaar (Haar wavelet), 
        DWTKernelpwl0, DWTKernelpwl2 (piecewise 
            linear wavelets with different number of van. moms.), 
        DWTKernel53 (Spline 5/3 wavelet, used for lossless 
            compression in JPEG2000),
        DWTKernel97 (CDF 9/7 wavelet, used for lossy compression in 
            JPEG2000),
        DWTKernelOrtho (Daubechies orthonormal wavelets with number 
            of vanishing moments dictated by a global variable)
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    for res in range(nres):
        f(x[0::2**res], symm, dual)
    reorganize_coefficients(x, nres, True)
            
def IDWTImpl(x, nres, f, symm=True, dual=False):
    """
    Compute the IDWT of x for a given number of resolutions, using 
    wavelet kernel f. The kernel is assumed to compute one level 
    of the IDWT in-place. IDWTImpl is responsible for reorganizing 
    the input so that the kernel can perform in-place calculation.
    The IDWT is computed along the first axis. x may have a second 
    axis, as is the case for sound with more than one channel. 
    The IDWT is then applied to each channel.
    
    x: The vector which we apply the IDWT to.
    nres: The number of stages
    f: The wavelet kernel to apply. Supported kernels are 
        IDWTKernelHaar (Haar wavelet), 
        IDWTKernelpwl0, DWTKernelpwl2 (piecewise 
            linear wavelets with different number of van. moms.), 
        IDWTKernel53 (Spline 5/3 wavelet, used for lossless 
            compression in JPEG2000),
        IDWTKernel97 (CDF 9/7 wavelet, used for lossy compression in 
            JPEG2000),
        IDWTKernelOrtho (Daubechies orthonormal wavelets with number 
            of vanishing moments dictated by a global variable)
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    reorganize_coefficients(x, nres, False)
    for res in range(nres - 1, -1, -1):
        f(x[0::2**res], symm, dual)
            
# Lifting steps
            
def liftingstepevensymm(lmbda, x, symm):
    """
    Apply an elementary symmetric lifting step of even type to x. 
    
    lmbda: The common value of the two filter coefficients
    x: The vector which we apply the lifting step to
    symm: Whether to apply symmetric extension to the input
    """
    if (not symm) and mod(len(x), 2)!=0:
        raise AssertionError()
    if symm:
        x[0] += 2*lmbda*x[1] # With symmetric extension
    else:
        x[0] += lmbda*(x[1]+x[-1])
    x[2:-1:2] += lmbda*(x[1:-2:2] + x[3::2])
    if mod(len(x), 2)==1 and symm:
        x[-1] += 2*lmbda*x[-2] # With symmetric extension
  
def liftingstepoddsymm(lmbda, x, symm):
    """
    Apply an elementary symmetric lifting step of odd type to x. 
    
    lmbda: The common value of the two filter coefficients
    x: The vector which we apply the lifting step to
    symm: Whether to apply symmetric extension to the input
    """
    if (not symm) and mod(len(x), 2)!=0:
        raise AssertionError()
    x[1:-1:2] += lmbda*(x[0:-2:2] + x[2::2])
    if mod(len(x), 2)==0:
        if symm:
            x[-1] += 2*lmbda*x[-2] # With symmetric extension
        else:
            x[-1] += lmbda*(x[0]+x[-2])

def liftingstepeven(lmbda1, lmbda2, x):
    """
    Apply an elementary non-symmetric lifting step of even type to x.
    
    lmbda1: The first filter coefficient
    lmbda2: The second filter coefficient
    x: The vector which we apply the lifting step to
    """
    if mod(len(x), 2)!=0:
        raise AssertionError()
    x[0] += lmbda1*x[1] + lmbda2*x[-1]
    x[2:-1:2] += lmbda1*x[3::2] + lmbda2*x[1:-2:2]
            
def liftingstepodd(lmbda1, lmbda2, x):
    """
    Apply an elementary non-symmetric lifting step of odd type to x.
    
    lmbda1: The first filter coefficient
    lmbda2: The second filter coefficient
    x: The vector which we apply the lifting step to
    """
    if mod(len(x), 2)!=0:
        raise AssertionError()
    x[1:-2:2] += lmbda1*x[2:-1:2] + lmbda2*x[0:-3:2]
    x[-1] += lmbda1*x[0] + lmbda2*x[-2]                                                

# The Haar wavelet

def DWTKernelHaar(x, symm, dual):
    """
    Apply the DWT kernel transformation for the Haar wavelet to x.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
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
         
def IDWTKernelHaar(x, symm, dual):
    """
    Apply the IDWT kernel transformation for the Haar wavelet to x.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
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
            
# Piecewise linear wavelets

def DWTKernelpwl0(x, symm, dual):
    """
    Apply the DWT kernel transformation for the 
    piecewise linear wavelet (i.e. 0 vanishing moments) to x.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    if dual:
        x /= sqrt(2)
        liftingstepevensymm(0.5, x, symm)
    else:
        x *= sqrt(2)
        liftingstepoddsymm(-0.5, x, symm)
        
def IDWTKernelpwl0(x, symm, dual):
    """
    Apply the IDWT kernel transformation for the 
    pieccewise linear wavelet (i.e. 0 vanishing moments) to x.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    if dual:
        x *= sqrt(2)
        liftingstepevensymm(-0.5, x, symm)
    else:
        x /= sqrt(2)
        liftingstepoddsymm(0.5, x, symm)


def DWTKernelpwl2(x, symm, dual):
    """
    Apply the DWT kernel transformation for the 
    alternative pieccewise linear wavelet (i.e. 2 van. moms.) to x.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    if dual:
        liftingstepevensymm(0.5, x, symm)
        liftingstepoddsymm(-0.25, x, symm)
        x /= sqrt(2)
    else:
        liftingstepoddsymm(-0.5, x, symm)
        liftingstepevensymm(0.25, x, symm)
        x *= sqrt(2)
    
def IDWTKernelpwl2(x, symm, dual):
    """
    Apply the IDWT kernel transformation for the 
    alternative pieccewise linear wavelet (i.e. 2 van. moms.) to x.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    if dual:
        x *= sqrt(2)
        liftingstepoddsymm(0.25, x, symm)
        liftingstepevensymm(-0.5, x, symm)
    else:
        x /= sqrt(2)
        liftingstepevensymm(-0.25, x, symm)
        liftingstepoddsymm(0.5, x, symm)       
                
        
# JPEG2000-related wavelet kernels
        
def DWTKernel53(x, symm, dual):
    """
    Apply the DWT kernel transformation for the 
    Spline 5/3 wavelet to x.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    if dual:
        x[0::2] *= 0.5
        x[1::2] *= 2
        liftingstepevensymm(0.125, x, symm)
        liftingstepoddsymm(-1, x, symm)
    else:
        x[0::2] *= 2
        x[1::2] *= 0.5
        liftingstepoddsymm(-0.125, x, symm)
        liftingstepevensymm(1, x, symm)
            
def IDWTKernel53(x, symm, dual):
    """
    Apply the IDWT kernel transformation for the 
    Spline 5/3 wavelet to x.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    if dual:
        liftingstepoddsymm(1, x, symm)
        liftingstepevensymm(-0.125, x, symm)     
        x[0::2] *= 2
        x[1::2] *= 0.5
    else:
        liftingstepevensymm(-1, x, symm)
        liftingstepoddsymm(0.125, x, symm)     
        x[0::2] *= 0.5
        x[1::2] *= 2


def DWTKernel97(x, symm, dual):
    """
    Apply the DWT kernel transformation for the CDF 9/7 wavelet to x.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    lambda1=-0.586134342059950
    lambda2=-0.668067171029734
    lambda3=0.070018009414994
    lambda4=1.200171016244178
    alpha=-1.149604398860250
    beta=-0.869864451624777
    if dual:
        x[0::2] /= alpha
        x[1::2] /= beta
        liftingstepevensymm(lambda4, x, symm)
        liftingstepoddsymm(lambda3, x, symm)
        liftingstepevensymm(lambda2, x, symm)
        liftingstepoddsymm(lambda1, x, symm)
    else:
        x[0::2] *= alpha
        x[1::2] *= beta
        liftingstepoddsymm(-lambda4, x, symm)
        liftingstepevensymm(-lambda3, x, symm)
        liftingstepoddsymm(-lambda2, x, symm)
        liftingstepevensymm(-lambda1, x, symm)
                
def IDWTKernel97(x, symm, dual):
    """
    Apply the IDWT kernel transformation for the CDF 9/7 wavelet to x
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    lambda1=-0.586134342059950
    lambda2=-0.668067171029734
    lambda3=0.070018009414994
    lambda4=1.200171016244178
    alpha=-1.149604398860250
    beta=-0.869864451624777
    if dual:
        liftingstepoddsymm(-lambda1, x, symm)
        liftingstepevensymm(-lambda2, x, symm)   
        liftingstepoddsymm(-lambda3, x, symm)
        liftingstepevensymm(-lambda4, x, symm)      
        x[0::2] *= alpha
        x[1::2] *= beta
    else:
        liftingstepevensymm(lambda1, x, symm)
        liftingstepoddsymm(lambda2, x, symm)   
        liftingstepevensymm(lambda3, x, symm)
        liftingstepoddsymm(lambda4, x, symm)      
        x[0::2] /= alpha
        x[1::2] /= beta
        
# Orthonormal wavelets
        
def DWTKernelOrtho( x, symm, dual):
    """
    Apply the DWT kernel transformation for orthonormal wavelets. 
    The number of vanishing moments is stored in global variables.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    global lambdas, alpha, beta
    global beta
    if dual:
        x[0::2] /= alpha
        x[1::2] /= beta
        for stepnr in range(lambdas.shape[0] - 1, 0, -2):
            liftingstepodd(lambdas[stepnr, 1], lambdas[stepnr, 0], x)
            liftingstepeven(lambdas[stepnr -1, 1], lambdas[stepnr - 1, 0], x)   
        if mod(lambdas.shape[0], 2)==1:
            liftingstepodd(lambdas[0, 1], lambdas[0, 0], x)
    else:
        x[0::2] *= alpha
        x[1::2] *= beta
        for stepnr in range(lambdas.shape[0] - 1, 0, -2):
            liftingstepeven(-lambdas[stepnr, 0], -lambdas[stepnr, 1], x)
            liftingstepodd(-lambdas[stepnr - 1, 0], -lambdas[stepnr - 1, 1], x)  
        if mod(lambdas.shape[0], 2)==1:
            liftingstepeven(-lambdas[0, 0], -lambdas[0, 1], x)
  
def IDWTKernelOrtho( x, symm, dual):
    """
    Apply the IDWT kernel transformation for orthonormal wavelets. 
    The number of vanishing moments is stored in global variables.
    
    x: The vector which we apply this kernel transformation to
    symm: Whether to apply symmetric extension to the input
    dual: Whether to apply the wavelet kernel or dual wavelet kernel.
    """
    global lambdas
    global alpha
    global beta
    if dual: 
        stepnr = 0
        if mod(lambdas.shape[0], 2) == 1: # Start with an odd step
            liftingstepodd(-lambdas[stepnr, 1], -lambdas[stepnr, 0], x)
            stepnr += 1
        while stepnr < lambdas.shape[0]:
            liftingstepeven(-lambdas[stepnr, 1], -lambdas[stepnr, 0], x)
            liftingstepodd(-lambdas[stepnr + 1, 1], -lambdas[stepnr + 1, 0], x)
            stepnr += 2
        x[0::2] *= alpha
        x[1::2] *= beta
    else:
        stepnr = 0
        if mod(lambdas.shape[0],2) == 1: # Start with an even step
            liftingstepeven(lambdas[stepnr, 0], lambdas[stepnr, 1], x)
            stepnr += 1
        while stepnr < lambdas.shape[0]:
            liftingstepodd(lambdas[stepnr, 0], lambdas[stepnr, 1], x)
            liftingstepeven(lambdas[stepnr + 1, 0], lambdas[stepnr + 1, 1], x)
            stepnr += 2
        x[0::2] /= alpha
        x[1::2] /=beta
            



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

def liftingfactortho(N):
    """
    Assume that len(h1)==len(h0), and that h0 and h1 are even length and as symmetric as possible, with h0 with a minimum possible overweight of filter coefficients to the left, h1 to the right
    This function computes lifting steps l1, l2,...,ln, and constants alpha, beta so that ln ... l2 l1 H =  diag(alpha,beta), and stores these as global variables
    This gives the following recipes for 
        Computing H: first multiply with diag(alpha,beta), then the inverses of the lifting steps in reverse order 
        Computing G: apply the lifting steps in the order given, finally multiply with diag(1/alpha,1/beta)
    ln is always odd, so that l1 is odd if and only if n is odd.
    All even lifting steps have only filter coefficients 0,1. All odd lifting steps have only filter coefficients -1,0
    """
    global lambdas, alpha, beta
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
    
def _test_kernel(f,invf,text):
    print text
    res = random.random(16)
    x = zeros(16)
    x[:] = res[:]
    DWTImpl(x,2,f)
    IDWTImpl(x,2,invf)
    diff = abs(x-res).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    res = random.random((16,2))
    x = zeros((16,2))
    x[:] = res[:]
    DWTImpl(x,2,f)
    IDWTImpl(x,2,invf)
    diff = abs(x-res).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
def _test_kernel_ortho():
    print 'Testing orthonormal wavelets'
    liftingfactortho(4)
    res = random.random(16) # only this assumes that N is even
    x = zeros(16)
    
    print 'Testing that the reverse inverts the forward transform'
    x[0:16] = res[0:16]
    DWTImpl(x, 2, DWTKernelOrtho)
    IDWTImpl(x, 2, IDWTKernelOrtho)
    diff = max(abs(x-res))
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing that the transform is orthogonal, i.e. that the transform and its dual are equal'
    x[0:16] = res[0:16]
    DWTImpl(x, 2, DWTKernelOrtho)
    DWTImpl(res, 2, DWTKernelOrtho, False, True)
    diff = max(abs(x-res))
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
def _test_dwt_different_sizes():
    print 'Testing the DWT on different input sizes'
    m = 4
    f = DWTKernel97
    invf = IDWTKernel97

    print 'Testing the DWT for greyscale image'
    img = random.random((32,32))
    img2 = zeros_like(img)
    img2[:] = img[:]
    DWT2Impl(img2, m, f)
    IDWT2Impl(img2, m, invf)
    diff = abs(img2-img).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing the DWT for RGB image'
    img = random.random((32, 32, 3))
    img2 = zeros_like(img)
    img2[:] = img[:]
    DWT2Impl(img2, m, f)
    IDWT2Impl(img2, m, invf)
    diff = abs(img2-img).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing the DWT for sound with one channel'
    sd = random.random(32)
    sd2 = zeros_like(sd)
    sd2[:] = sd[:]
    DWTImpl(sd2, m, f)
    IDWTImpl(sd2, m, invf)
    diff = abs(sd2-sd).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing the DWT for sound with two channels'
    sd = random.random((32,2))
    sd2 = zeros_like(sd)
    sd2[:] = sd[:]
    DWTImpl(sd2, m, f)
    IDWTImpl(sd2, m, invf)
    diff = abs(sd2-sd).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
def _test_orthogonality():
    print 'Testing that the wavelet and the dual wavelet are equal for orthonormal wavelets'
    liftingfactortho(4)
    x0 = random.random(32)
    
    print 'Testing that the IDWT inverts the DWT'
    x = x0.copy()
    DWTImpl(x, 2, DWTKernelOrtho, 0, 0)
    IDWTImpl(x, 2, IDWTKernelOrtho, 0, 0);
    diff = abs(x-x0).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Apply the transpose, to see that the transpose equals the inverse'
    x = x0.copy()
    DWTImpl(x, 2, DWTKernelOrtho, 0, 0)
    IDWTImpl(x, 2, IDWTKernelOrtho, 0, 1)
    diff = abs(x-x0).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff

    print 'To see this at the level of kernel transformations'
    x = x0.copy()
    DWTKernelOrtho(x, 0, 0)
    IDWTKernelOrtho(x, 0, 1)
    diff = abs(x-x0).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff

    print 'See that the wavelet transform equals the dual wavelet transform'
    x = x0.copy()
    DWTImpl(x, 2, DWTKernelOrtho, 0, 1)
    DWTImpl(x0, 2, DWTKernelOrtho, 0, 0)
    diff = abs(x-x0).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff



if __name__=='__main__':
    _test_dwt_different_sizes()
    _test_kernel_ortho()
    _test_orthogonality()
    _test_kernel( DWTKernel97, IDWTKernel97, 'Testing CDF 9/7 wavelet')
    _test_kernel( DWTKernel53, IDWTKernel53, 'Testing Spline 5/3 wavelet')
    _test_kernel( DWTKernelpwl0, IDWTKernelpwl0, 'Testing piecewise linear wavelet')
    _test_kernel( DWTKernelpwl2, IDWTKernelpwl2, 'Testing alternative piecewise linear wavelet')
    _test_kernel( DWTKernelHaar, IDWTKernelHaar, 'Testing Haar wavelet')