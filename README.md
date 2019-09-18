# Linear Algebra, Signal Processing, and Wavelets - A Unified Approach

This repository contains all the source code related to the book *"Linear
Algebra, Signal Processing, and Wavelets - A Unified Approach"* (both MatLab and
Python version), by Øyvind Ryan, Springer, 2019. 
It contains open source implementations of the fast Fourier transform (FFT), 
discrete cosine transform (DCT) and the discrete wavelet transform (DWT).
The implementations are educational in nature. In particular, although efficient code is attempted, efficiency is sacrificed some places in favor of clarity of code. 
In particular we highlight the advantages of the implemented DWT.

* While Spline wavelets use a filter-based DWT implemention, orthonormal wavelets based on the construction by Daubechies, as well as some least dissimilar symmetric wavelets, use a 
  lifting-based implementation. Lifting-based implementations roughly reduce the number of operations by a factor of 2.
  A change in the API where the caller can indicate either a lifting- or filter-based implementation is under preparation. This will be useful for testing purposes, as it will let the caller compare the two approaches.
* The API currently supports wavelets preserving vanishing moments on the interval, for Spline wavelets and for orthonormal wavelets. A paper is submitted which explains the generalization to other wavelets, and a revision of the source code 
  is under way which will capture the results therein.
  To the best of our knowledge this is the only openly available software implementation which computes boundary wavelets on the fly..  It also works for general input sizes, and it turns out that the boundary wavelets themselves depend heavily on the input size.)

## Installation
### Python 
Add the `python` directory to your python path environment variable e.g. on
Unix systems you can type 
```
export PYTHONPATH=$PYTHONPATH:/path/to/applinalgcode/python
```

### MatLab
If you need to run examples from the book  
you need to add the folders `matlab`, `images` and `sounds` to your MatLab path. It can be an
advantage to add this automatically on startup. In this case add it to your 
[startup.m file](https://ch.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html). We have created the file `install.m` which will try to do this automatically 
for you.

## User guide to the discrete wavelet transform implementation

The software has two main functions `dwt_impl` and `idwt_impl` which implement
a DWT and a IDWT, respectively. These functions can be called with a variety of
arguments, specifying the type of wavelet transform and its properties. To
compute an `m`-level DWT on a vector `x`, using a wavelet with name `wname`,
boundary handling mode `bd_mode`, and prefiltering mode `prefilter_mode`,
simply write
```
y = dwt_impl(x, wname, m, bd_mode, prefilter_mode);
```
Here `wname` will be one of the following supported wavelets:
*  `'cdf97'` : CDF 9/7 wavelet with $N = \tilde{N} = 4$ vanishing moments.
*  `'spline53'` : Spline 5/3 wavelet with $N = \tilde{N} = 2$ vanishing moments.
*  `'splineX.X'` : Spline wavelet with `X` number of vanishing moments for the wavelet and dual wavelet,
*  `'pwl0'` : Piecewise linear wavelets with 0 vanishing moments,
*  `'pwl2'` : Piecewise linear wavelets with 2 vanishing moments,
*  `'Haar'` : The Haar wavelet,
*  `'dbX'` : Daubechies orthonormal wavelet with `X` vanishing moments,
*  `'symX'` : Symlets. A close to symmetric, orthonormal (Daubechies) wavelet with `X` vanishing moments.

Likewise `bd_mode` can take the values
* `'per'` : Periodic extension,
* `'symm'` : Symmetric extension (not for orthonormal wavelets),
* `'none'` : No boundary handling, in the sense that the input is zero-padded. 
* `'bd'` : Boundary handling as described in XXXX (Only available in MatLab version) 

Note that, when the `'bd'` mode is invoked, computations in the `'none'` mode are also performed by the software, in order to compute the tail handling components. 
`prefilter_mode` can take the values
* `'none'` : No prefiltering (default),
* `'bd_pre'` : Boundary wavelets with preconditioning as described in XXXX (Only available in MatLab version) 

Not all combinations of these arguments make sense. For instance it is not possible to apply a symmetric boundary extension to an orthonormal
wavelet. In such cases the functions halt, issuing an error. 

`dwt_impl` also accepts the following arguments.
* `dims` : The number of dimensions to apply the DWT to. If the input is two-dimensional, this enables the caller to specify whether a two dimensional DWT should be applied, or a one dimenionsal DWT vectorized on the second axis. 
* `dual` :  Whether the dual transform should be applied or not.
* `transpose` : Whether the transform shoudl be transposed.
* `data_layout` : How data should be assembled. Possible values are `resolution` (lowest resolution first, default), and `time` (sort according to time).


## Internal functions and efficient computations 
The `dwt_impl` and `idwt_impl` functions compute the
filter coefficients and tail handling components on the fly each time the
functions are invoked. This allows the software to compute different boundary
functions for different input sizes and makes the functions user friendly. At
the same time, this increases the computational time as each call recomputes coefficients. To allow for using precomputed boundary functions, the software has the functions `dwt1_impl_internal` and
`idwt1_impl_internal` which will (given the right input) compute a one dimensional DWT and IDWT, respectively, using precomputed boundary functions.  Similar functions exist for two and three
dimensions. 

To use these internal functions we need to first compute the filter coefficients 
manually. The complete code can be as follows.
```
[wav_props, dual_wav_props] = find_wav_props(m, wname, bd_mode, size(x,1));
[f, prefilter] = find_kernel(wav_props, dual_wav_props, forward, dual, ...
                             transpose, prefilter_mode);

offsets = [wav_propsx.offset_L, wav_propsx.offset_R];

y = dwt1_impl_internal(x, f, m, bd_mode, prefilter, offsets);
```
Let us go through the different pieces of this code. 
* `wav_props` and `dual_wav_props`, as returned by `find_wav_props`, are structures (objects in Python) which contain all filter coefficients and matrices related to the boundary handling. 
* The `find_kernel` functions return two function handles:
  1. `f(x, bd_mode)` preforms a one level DWT on `x` with boundary mode `bd_mode`, 
  2. `prefilter(x, forward)` filters `x`, where `forward` is either 0 (postfiltering) or 1 (prefiltering). 
* The `offset` parameter is only necessary when using boundary wavelets. It tells `dwt1_impl_internal` that there are $N-K_L$ and $N-K_R$ extra wavelets at each boundary.

Since the internal functions avoid precomputation, their execution times 
should be comparable with those of convolution, since the DWT/IDWT can be 
expressed in terms of this (and possibly lifting).

Following the code above, it is possible to experiment with custom made kernels - simply implement your own kernel function (taking the same arguments as the ones defined in the software), and use it as input to `dwt1_impl_internal`.




