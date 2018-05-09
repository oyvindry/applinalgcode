function x=DWTImpl_internal(x, m, dwt_kernel, bd_mode, prefilter, wav_props, data_layout)
    inds = 1:size(x,1);
    x = prefilter(x, 1);
    for res=0:(m - 1)
        x(inds, :) = dwt_kernel(x(inds, :), bd_mode);
        inds = inds((wav_props.offset_L+1):2:(end-wav_props.offset_R));
    end
    x(inds, :) = prefilter(x(inds, :),0);
    x = reorganize_coeffs_forward(x, m, wav_props, data_layout);
end

function y=reorganize_coeffs_forward(x, m, wav_props, data_layout)
    y = x;
    if strcmpi(data_layout, 'resolution')
        N = size(x,1);
        inds = 1:N;
        endy = N;
        for res=1:m
            xindices = [inds(1:wav_props.offset_L) inds((wav_props.offset_L + 2):2:(end-wav_props.offset_R)) inds((end-wav_props.offset_R+1):end)]; % psi-indices
            y((endy-length(xindices)+1):endy,:) = x(xindices,:);
            endy = endy-length(xindices);
            inds = inds((wav_props.offset_L+1):2:(end-wav_props.offset_R)); 
        end
        y(1:endy, :) = x(inds, :);
    end
end