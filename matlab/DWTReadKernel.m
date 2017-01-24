function x = DWTReadKernel(x, nres, f, mode, dualarg)
    % DWTReadKernel computes the DWT of the vector x for a given number of
    % resolutions, using a given wavelet kernel. The kernel is assumed to
    % compute one level of the DWT, using the order dictated by in-place
    % computation. DWTReadKernel is responsible for reorganizing the output so
    % that the low resolution coefficients comes first, as required by the DWT.
    % The DWT is computed along the first dimension.  x may have a second
    % dimension, as is the case for sound with more than one channel. The DWT is
    % then applied to each channel.
    %   
    % x: The vector which we apply the DWT to.
    % nres: The number of stages
    % f: The wavelet kernel to apply. Supported kernels are 
    %     DWTKernelHaar (Haar wavelet), 
    %     DWTKernelpwl0, DWTKernelpwl2 (piecewise linear wavelets with 
    %                                   different number of vanishing moments), 
    %     DWTKernel53 (Spline 5/3 wavelet, 
    %                 used for lossless compression in JPEG2000),
    %     DWTKernel97 (CDF 9/7 wavelet, used for lossy compression in JPEG2000),
    %     DWTKernelOrtho (Orthonormal wavelets with a given number of vanishing 
    %                     moments)
    % mode: Boundary extension mode
    % dual: Whether to apply the kernel for the dual wavelet rather than the 
    %       wavelet itself.
    
    N = size(x, 1);
    for res=0:(nres - 1)
        x(1:2^res:N, :) = f(x(1:2^res:N, :), mode, dualarg);
    end
    x = reorganize_coefficients(x, nres, 1);
end
