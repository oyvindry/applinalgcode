from math import *
from numpy import *
from sound import *
from images import *
from scipy.fftpack import dct

#bitreverse

def bitreverse(x):
    """
    bitreverse(x) alters x in bit-reversed order.
    """
    N = len(x)
    j = 0
    for i in range(0, N/2, 2):
        if (j > i):         
            temp = x[j]; x[j] = x[i]; x[i] = temp
            temp = x[N - j - 1]; x[N - j - 1] = x[N - i - 1]; x[N - i - 1] = temp
        temp = x[i+1]; x[i+1] = x[j + N/2]; x[j + N/2] = temp
        m = N/4
        while (m >= 1 and j >= m):
            j -= m
            m /= 2
        j += m
    
    
def bitreversearr(x):
    N, n = shape(x)
    temp=zeros(n).astype(complex)
    j = 0
    for i in range(0, N/2, 2):
        if (j > i):         
            temp[:] = x[j]; x[j, :] = x[i]; x[i, :] = temp
            temp[:] = x[N - j - 1]; x[N - j - 1, :] = x[N - i - 1]; x[N - i - 1, :] = temp
        temp[:] = x[i+1]; x[i+1, :] = x[j + N/2]; x[j + N/2, :] = temp    
        m = N/4
        while (m >= 1 and j >= m):
            j -= m
            m /= 2
        j += m

# DCT code

def DCTImpl(x):
    """ 
    Compute the DCT of the vector x
    
    x: a vector
    """ 
    N = len(x)
    if N > 1:
        x1 = concatenate([x[0::2], x[-1:0:-2]]).astype(complex)
        FFTImpl(x1, FFTKernelStandard)
        cosvec = cos(pi*arange(float(N))/(2*N))
        sinvec = sin(pi*arange(float(N))/(2*N))
        if ndim(x) == 1:
            x[:] = cosvec*real(x1) + sinvec*imag(x1)
        else: 
            for s2 in xrange(shape(x)[1]):
                x[:, s2] = cosvec*real(x1[:, s2]) \
                + sinvec*imag(x1[:, s2])
        x[0] *= sqrt(1/float(N))
        x[1:] *= sqrt(2/float(N))
        
def IDCTImpl(y):
    """ 
    Compute the IDCT of the vector y
    
    y: a vector
    """
    N = len(y)
    if N > 1:
        y[0] /= sqrt(1/float(N))
        y[1:] /= sqrt(2/float(N))
        Q = exp(-pi*1j*arange(float(N))/(2*N))
        y1 = zeros_like(y).astype(complex)
        y1[0] = y[0]/Q[0]
        if ndim(y) == 1:
            y1[1:] = (y[1:] - 1j*y[-1:0:-1])/Q[1:]
        else:
            for s2 in xrange(shape(y)[1]):
                y1[1:, s2] = (y[1:, s2] - 1j*y[-1:0:-1, s2])/Q[1:]
        FFTImpl(y1, FFTKernelStandard, 0)
        y[0::2] = real(y1[0:(N/2)])
        y[1::2] = real(y1[-1:(N/2-1):-1])

def DFTImpl(x, forward=True):
    """
    Compute the DFT of the vector x using standard matrix 
    multiplication. To avoid out of memory situations, we do not 
    allocate the entire DFT matrix, only one row of it at a time. 
    Note that this function differs from the FFT in that it includes 
    the normalizing factor 1/sqrt(N). The DFT is computed along axis 
    0. If there is another axis, the DFT is computed for each element
    in this as well. 
    
    x: a vector
    forward: Whether or not this is forward (i.e. DFT) 
    or reverse (i.e. IDFT)
    """
    y = zeros_like(x).astype(complex)
    N = len(x)
    sign = -(2*forward - 1) 
    if ndim(x) == 1:
        for n in xrange(N):
            D = exp(sign*2*pi*n*1j*arange(float(N))/N)
            y[n] = dot(D, x)
    else:
        for n in range(N):
            D = exp(sign*2*pi*n*1j*arange(float(N))/N)
            for s2 in xrange(shape(x)[1]):
                y[n,s2] = dot(D,x[:, s2])
    if sign == 1:
        y /= float(N)
    return y

def FFTImpl(x, FFTKernel, forward = True):
    """
    Compute the FFT or IFFT of the vector x. Note that this function 
    differs from the DFT in that the normalizing factor 1/sqrt(N) is 
    not included. The FFT is computed along axis 0. If there is 
    another axis, the FFT is computed for each element in this as 
    well. This function calls a kernel for computing the FFT. The 
    kernel assumes that the input has been bit-reversed, and contains
    only one axis. This function is where the actual bit reversal and 
    the splitting of the axes take place.
    
    x: a vector
    FFTKernel: can be any of FFTKernelStandard, FFTKernelNonrec, and 
    FFTKernelSplitradix. The kernel assumes that the input has been 
    bit-reversed, and contains only one axis. 
    forward: Whether the FFT or the IFFT is applied
    """
    if ndim(x) == 1:
        bitreverse(x)
        FFTKernel(x, forward)
    else:
        bitreversearr(x)
        for s2 in xrange(shape(x)[1]):
            FFTKernel(x[:, s2], forward)
    if not forward:
        x /= len(x)

def FFTKernelStandard(x, forward):
    """
    Compute the FFT of x, using a standard FFT algorithm.
    
    x: a bit-reversed version of the input. Should have only one axis
    forward: Whether the FFT or the IFFT is applied
    """
    N = len(x)
    sign = -1
    if not forward:
        sign = 1
    if N > 1:
        xe, xo = x[0:(N/2)], x[(N/2):]
        FFTKernelStandard(xe, forward)
        FFTKernelStandard(xo, forward)
        D = exp(sign*2*pi*1j*arange(float(N/2))/N) 
        xo *= D 
        x[:] = concatenate([xe + xo, xe - xo]) 
                
def FFTKernelNonrec(x, forward):
    """
    Compute the FFT of x, using a non-recursive FFT algorithm. 
    
    x: a bit-reversed version of the input. Should have only one axis
    forward: Whether the FFT or the IFFT is applied
    """
    N = len(x)
    sign = -1
    if not forward:
        sign = 1
    D = exp(sign*2*pi*1j*arange(float(N/2))/N)
    nextN = 1
    while nextN < N:
        k = 0
        while k < N:
            xe, xo = x[k:(k + nextN)], x[(k + nextN):(k + 2*nextN)]
            xo *= D[0::(N/(2*nextN))]
            x[k:(k+2*nextN)] = concatenate([xe + xo, xe - xo])
            k += 2*nextN
        nextN *= 2     
    
def FFTKernelSplitradix(x, forward):
    """
    Compute the FFT of x, using the split-radix FFT algorithm. 
    
    x: a bit-reversed version of the input. Should have only one axis
    forward: Whether the FFT or the IFFT is applied
    """
    N = len(x)
    sign = -1
    if not forward:
        sign = 1
    if N == 2:
        x[:] = [x[0] + x[1], x[0] - x[1]]
    elif N > 2:
        xe, xo1, xo2 = x[0:(N/2)], x[(N/2):(3*N/4)], x[(3*N/4):N]
        FFTKernelSplitradix(xe, forward)
        FFTKernelSplitradix(xo1, forward)
        FFTKernelSplitradix(xo2, forward)
        G = exp(sign*2*pi*1j*arange(float(N/4))/N)
        H = G*exp(sign*2*pi*1j*arange(float(N/4))/(N/2))
        xo1 *= G
        xo2 *= H
        xo = concatenate( [xo1 + xo2, -sign*1j*(xo2 - xo1)] )
        x[:] = concatenate([xe + xo, xe - xo])   
    
def tensor_impl(X, S1, S2):
    M, N = shape(X)[0:2]
    for n in range(N):
        S1(X[:,n])
    for m in range(M):
        S2(X[m, :])
    
# We can use classes to make a more efficient FFT implementation
    
def _testfft():   
    x1 = random.random(32).astype(complex)
    x2 = random.random((32,3)).astype(complex)
    x2copy = zeros_like(x2); x2copy[:] = x2[:]
    
    print 'Testing DFTImpl'
    x = fft.fft(x1, axis=0)
    x0 = DFTImpl(x1)
    diff = abs(x0-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x = fft.fft(x2, axis=0)
    x0 = DFTImpl(x2)
    diff = abs(x0-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing FFTKernelStandard'
    x = fft.fft(x1, axis=0)
    x1copy = x1.copy()
    FFTImpl(x1copy, FFTKernelStandard) 
    diff = abs(x1copy-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x = fft.fft(x2, axis=0)
    x2copy = x2.copy()
    FFTImpl(x2copy, FFTKernelStandard) 
    diff = abs(x2copy-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing FFTKernelNonrec'
    x = fft.fft(x1, axis=0)
    x1copy = x1.copy()
    FFTImpl(x1copy, FFTKernelNonrec) 
    diff = abs(x1copy-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x = fft.fft(x2, axis=0)
    x2copy = x2.copy()
    FFTImpl(x2copy, FFTKernelNonrec) 
    diff = abs(x2copy-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing FFTKernelSplitradix'
    x = fft.fft(x1, axis=0)
    x1copy = x1.copy()
    FFTImpl(x1copy, FFTKernelSplitradix) 
    diff = abs(x1copy-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x = fft.fft(x2, axis=0)
    x2copy = x2.copy()
    FFTImpl(x2copy, FFTKernelSplitradix) 
    diff = abs(x2copy-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing that FFTKernelStandard(reverse) inverts FFTKernelStandard(forward)' 
    x1copy = x1.copy()
    FFTImpl(x1copy, FFTKernelStandard) 
    FFTImpl(x1copy, FFTKernelStandard, 0) 
    diff = abs(x1copy-x1).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x2copy = x2.copy()
    FFTImpl(x2copy, FFTKernelStandard) 
    FFTImpl(x2copy, FFTKernelStandard, 0) 
    diff = abs(x2copy-x2).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing that FFTKernelNonrec(reverse) inverts FFTKernelNonrec(forward)' 
    x1copy = x1.copy()
    FFTImpl(x1copy, FFTKernelNonrec) 
    FFTImpl(x1copy, FFTKernelNonrec, 0) 
    diff = abs(x1copy-x1).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x2copy = x2.copy()
    FFTImpl(x2copy, FFTKernelNonrec) 
    FFTImpl(x2copy, FFTKernelNonrec, 0) 
    diff = abs(x2copy-x2).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing that FFTKernelSplitradix(reverse) inverts FFTKernelSplitradix(forward)' 
    x1copy = x1.copy()
    FFTImpl(x1copy, FFTKernelSplitradix) 
    FFTImpl(x1copy, FFTKernelSplitradix, 0) 
    diff = abs(x1copy-x1).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x2copy = x2.copy()
    FFTImpl(x2copy, FFTKernelSplitradix) 
    FFTImpl(x2copy, FFTKernelSplitradix, 0) 
    diff = abs(x2copy-x2).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    print 'Testing that DFTImpl(reverse) inverts DFTImpl(forward)'
    x0 = DFTImpl(DFTImpl(x1), 0)
    diff = abs(x0-x1).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x0 = DFTImpl(DFTImpl(x2), 0)
    diff = abs(x0-x2).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff

def _testdct():
    x1 = random.random(32)
    x2 = random.random((32,3))
    x2copy = zeros_like(x2); x2copy[:] = x2[:]
    
    print 'Testing DCTImpl'
    x = dct(x1, norm='ortho', axis=0)
    x1copy = x1.copy()
    DCTImpl(x1copy) 
    diff = abs(x1copy-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x = dct(x2, norm='ortho', axis=0)
    x2copy = x2.copy()
    DCTImpl(x2copy) 
    diff = abs(x2copy-x).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    
    print 'Testing that IDCTImpl inverts DCTImpl'
    x1copy = x1.copy()
    DCTImpl(x1copy) 
    IDCTImpl(x1copy) 
    diff = abs(x1copy-x1).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
    x2copy = x2.copy()
    DCTImpl(x2copy) 
    IDCTImpl(x2copy) 
    diff = abs(x2copy-x2).max()
    assert diff < 1E-13, 'bug, diff=%s' % diff
    
if __name__=='__main__':
    _testfft()
    _testdct()
    
# For classes: Improvements
# 1. Stores tables for cosine/sine values    