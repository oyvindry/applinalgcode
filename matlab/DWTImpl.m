function x = DWTImpl(x, nres, wave_name, mode, dualarg)
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
    %            'dbX'   - Dauberchies orthnormal wavelet with X vanising
    %                      moments
    %            'symX'  - Symmlets: A close to symmetric, orthnormal wavelet 
    %                      with X vanising moments
    %
    % mode: Boundary extension mode
    % dual: Whether to apply the kernel for the dual wavelet rather than the 
    %       wavelet itself.
    % 
    if (~exist('mode')) mode = 1; end
    if (~exist('dualarg')) dualarg  = 0; end
    
    f = findDWTKernel(wave_name);
    x = DWTReadKernel(x, nres, f, mode, dualarg);
end
