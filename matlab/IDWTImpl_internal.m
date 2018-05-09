function x=IDWTImpl_internal(x, m, f, bd_mode, prefilter, wav_props, data_layout)
    [x, resstart, resend] = reorganize_coeffs_reverse(x, m, wav_props, data_layout);
    x(resstart(m+1):2^m:resend(m+1)) = prefilter(x(resstart(m+1):2^m:resend(m+1)), 1);
    for res = (m - 1):(-1):0
        inds = resstart(res+1):2^res:resend(res+1);
        x(inds, :) = f(x(inds, :), bd_mode);
    end
    x = prefilter(x, 0);
end

function [y, resstart, resend]=reorganize_coeffs_reverse(x, m, wav_props, data_layout)
    N = size(x,1);
    inds = 1:N;
    y = x;
    resstart = 1:(m+1); resend = 1:(m+1);
    resstart(1) = inds(1);
    resend(1) = inds(end);
    if strcmpi(data_layout, 'time')
        for res=1:m
            inds = inds((wav_props.offset_L+1):2:(end-wav_props.offset_R)); 
            resstart(res+1) = inds(1);
            resend(res+1) = inds(end);
        end
    end
    if strcmpi(data_layout, 'resolution')
        endy = N;
        for res=1:m
            xindices = [inds(1:wav_props.offset_L) inds((wav_props.offset_L + 2):2:(end-wav_props.offset_R)) inds((end-wav_props.offset_R+1):end)]; % psi-indices
            resstart(res+1) = inds(wav_props.offset_L+1);
            resend(res+1)   = inds(end-wav_props.offset_R);
            y(xindices,:) = x((endy-length(xindices)+1):endy,:);
            endy = endy-length(xindices);
            inds = inds((wav_props.offset_L+1):2:(end-wav_props.offset_R));
        end
        y(inds,:) = x(1:endy,:);
    end
end