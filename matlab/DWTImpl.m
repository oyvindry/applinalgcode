function x = DWTImpl(x, nres, wave_name, mode, dual)
    % Compute the discrete wavelet transform of the matrix `x` along the first
    % dimension.
    %    
    % x:         Matrix whose DWT will be computed along the first dimension        
    % nres:      Number of wavelet decompositions
    % wave_name: Name of the wavelet. Possible arguments are:
    %            'cdf97' - CDF 9/7 wavelet
    %            'cdf53' - Spline 5/3 wavelet  
    %            'pwl0'  - Piecewise linear wavelets with 0 vanishing moments
    %            'pwl2'  - Piecewise linear wavelets with 2 vanishing moments
    %            'Haar'  - The Haar wavelet
    %            'dbX'   - Dauberchies orthnormal wavelet with X vanishing
    %                      moments
    %            'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
    %                      with X vanishing moments
    %
    % mode: Boundary extension mode
    % dualarg: Whether to apply the kernel for the dual wavelet rather than the 
    %       wavelet itself.
    % 
    if (~exist('mode')) mode = 1; end
    if (~exist('dual')) dual  = 0; end
    
    f = findDWTKernel(wave_name);
    for res=0:(nres - 1)
        x(1:2^res:end, :) = f(x(1:2^res:end, :), mode, dual);
    end
    x = reorganize_coefficients(x, nres, 1);
end
