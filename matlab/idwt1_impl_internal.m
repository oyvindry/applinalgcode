function x=idwt1_impl_internal(x, f, m, bd_mode, prefilter, offsets, data_layout)
    % Compute a 1D IDWT using a precomputed kernel. The kernel may be the default library kernel obtained by calling find_kernel, 
    % or may be used-defined.
    %
    % x:         Matrix whose IDWT will be computed along the first dimension(s). 
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
    
    [x, resstart, resend] = reorganize_coeffs_reverse(x, m, offsets, data_layout);
    x(resstart(m+1):2^m:resend(m+1)) = prefilter(x(resstart(m+1):2^m:resend(m+1)), 1);
    for res = (m - 1):(-1):0
        inds = resstart(res+1):2^res:resend(res+1);
        x(inds, :) = f(x(inds, :), bd_mode);
    end
    x = prefilter(x, 0);
end

function [y, resstart, resend]=reorganize_coeffs_reverse(x, m, offsets, data_layout)
    N = size(x,1);
    inds = 1:N;
    y = x;
    resstart = 1:(m+1); resend = 1:(m+1);
    resstart(1) = inds(1);
    resend(1) = inds(end);
    if strcmpi(data_layout, 'time')
        for res=1:m
            inds = inds((offsets(1,1)+1):2:(end-offsets(1,2)));
            resstart(res+1) = inds(1);
            resend(res+1) = inds(end);
        end
    end
    if strcmpi(data_layout, 'resolution')
        endy = N;
        for res=1:m
            xindices = [inds(1:offsets(1,1)) inds((offsets(1,1) + 2):2:(end-offsets(1,2))) inds((end-offsets(1,2)+1):end)]; % psi-indices
            resstart(res+1) = inds(offsets(1,1)+1);
            resend(res+1)   = inds(end-offsets(1,2));
            y(xindices,:) = x((endy-length(xindices)+1):endy,:);
            endy = endy-length(xindices);
            inds = inds((offsets(1,1)+1):2:(end-offsets(1,2)));
        end
        y(inds,:) = x(1:endy,:);
    end
end