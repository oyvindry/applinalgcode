function x=dwt1_impl_internal(x, f, m, bd_mode, prefilter, offsets, data_layout)
    % Compute a 1D DWT using a precomputed kernel. The kernel may be the default library kernel obtained by calling find_kernel, 
    % or may be used-defined.
    %
    % x:         Matrix whose DWT will be computed along the first dimension(s). 
    % f:         kernel function     
    % m:         Number of resolutions. Default is 1
    % bd_mode:   Boundary extension mode. Possible modes are. 
    %            'per'    - Periodic extension
    %            'symm'   - Symmetric extension (default)
    %            'none'   - Take no extra action at the boundaries
    %            'bd'     - Boundary wavelets
    % prefilter: function which computes prefiltering. The default is no prefiltering.
    % offsets:   offsets at the beginning and the end as used by boundary wavelets. Default: zeros.
    % data_layout: How data should be assembled. Possible modes are:
    %            'resolution': Lowest resolution first (default)
    %            'time': Sort according to time
    
    if (~exist('m','var')) m = 1; end
    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('prefilter','var')) prefilter = @(x, forward) x; ; end
    if (~exist('offsets','var')) offsets = zeros(1,2); end
    if (~exist('data_layout','var')) data_layout = 'resolution'; end
    
    inds = 1:size(x,1);
    x = prefilter(x, 1);
    for res=0:(m - 1)
        x(inds, :) = f(x(inds, :), bd_mode);
        inds = inds((offsets(1,1)+1):2:(end-offsets(1,2)));
    end
    x(inds, :) = prefilter(x(inds, :),0);
    x = reorganize_coeffs_forward(x, m, offsets, data_layout);
end

function y=reorganize_coeffs_forward(x, m, offsets, data_layout)
    y = x;
    if strcmpi(data_layout, 'resolution')
        N = size(x,1);
        inds = 1:N;
        endy = N;
        for res=1:m
            xindices = [inds(1:offsets(1,1)) inds((offsets(1,1) + 2):2:(end-offsets(1,2))) inds((end-offsets(1,2)+1):end)]; % psi-indices
            y((endy-length(xindices)+1):endy,:) = x(xindices,:);
            endy = endy-length(xindices);
            inds = inds((offsets(1,1)+1):2:(end-offsets(1,2))); 
        end
        y(1:endy, :) = x(inds, :);
    end
end