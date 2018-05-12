function x=IDWT3Impl_internal(x, m, fx, fy, fz, bd_mode, prefilterx, prefiltery, prefilterz, wav_propsx, wav_propsy, wav_propsz, data_layout)
    [x, resstart, resend] = reorganize_coeffs3_reverse(x, m, wav_propsx, wav_propsy, wav_propsz, data_layout);
    
    lastdim = 1;
    if length(size(x)) == 4
        lastdim = size(x, 4);
    end 
    
    % postconditioning
    indsx = resstart(1,m+1):2^m:resend(1,m+1); indsy = resstart(2,m+1):2^m:resend(2,m+1); indsz = resstart(3,m+1):2^m:resend(3,m+1);
    x=tensor3_kernel(x, indsx, indsy, indsz, @(x,bd_mode) prefilterx(x, 1), @(x,bd_mode) prefiltery(x, 1), @(x,bd_mode) prefilterz(x, 1), lastdim, bd_mode);

    for res = (m - 1):(-1):0
        indsx = resstart(1,res+1):2^res:resend(1,res+1); 
        indsy = resstart(2,res+1):2^res:resend(2,res+1);
        indsz = resstart(3,res+1):2^res:resend(3,res+1);
        x = tensor3_kernel(x, indsx, indsy, indsz, fx, fy, fz, lastdim, bd_mode);
    end
    
    % preconditioning
    indsx = resstart(1,1):resend(1,1); indsy = resstart(2,1):resend(2,1); indsz = resstart(3,1):resend(3,1);
    x = tensor3_kernel(x, indsx, indsy, indsz, @(x,bd_mode) prefilterx(x, 0), @(x,bd_mode) prefiltery(x, 0), @(x,bd_mode) prefilterz(x, 0), lastdim, bd_mode);
end

function [sig_out, resstart, resend]=reorganize_coeffs3_reverse(sig_in, m, wav_propsx, wav_propsy, wav_propsz, data_layout)
    indsx = 1:size(sig_in,1); indsy = 1:size(sig_in,2); indsz = 1:size(sig_in,3);
    sig_out = sig_in;
    resstart = [1:(m+1); 1:(m+1); 1:(m+1)]; resend = [1:(m+1); 1:(m+1); 1:(m+1)];
    resstart(1,1) = indsx(1); resend(1,1) = indsx(end);
    resstart(2,1) = indsy(1); resend(2,1) = indsy(end);
    resstart(3,1) = indsz(1); resend(3,1) = indsz(end);
    if strcmpi(data_layout, 'time')
        for res=1:m
            indsx = indsx((wav_propsx.offset_L+1):2:(end-wav_propsx.offset_R));
            indsy = indsx((wav_propsy.offset_L+1):2:(end-wav_propsy.offset_R));
            indsz = indsz((wav_propsz.offset_L+1):2:(end-wav_propsz.offset_R));
            resstart(1,res+1) = indsx(1); resend(1,res+1) = indsx(end);
            resstart(2,res+1) = indsy(1); resend(2,res+1) = indsy(end);
            resstart(3,res+1) = indsz(1); resend(3,res+1) = indsz(end);
        end
    end
    if strcmpi(data_layout, 'resolution')
        endx = size(sig_in,1); endy = size(sig_in,2); endz = size(sig_in,3);
        for res=1:m
            psiinds_x = [indsx(1:wav_propsx.offset_L) indsx((wav_propsx.offset_L + 2):2:(end-wav_propsx.offset_R)) indsx((end-wav_propsx.offset_R+1):end)]; % psi-indices
            psiinds_y = [indsy(1:wav_propsy.offset_L) indsy((wav_propsy.offset_L + 2):2:(end-wav_propsy.offset_R)) indsy((end-wav_propsy.offset_R+1):end)];
            psiinds_z = [indsz(1:wav_propsz.offset_L) indsz((wav_propsz.offset_L + 2):2:(end-wav_propsz.offset_R)) indsz((end-wav_propsz.offset_R+1):end)];
            phiinds_x = indsx((wav_propsx.offset_L + 1):2:(end-wav_propsx.offset_R));
            phiinds_y = indsx((wav_propsy.offset_L + 1):2:(end-wav_propsy.offset_R));
            
            resstart(1,res+1) = indsx(wav_propsx.offset_L+1); resend(1,res+1)   = indsx(end-wav_propsx.offset_R);
            resstart(2,res+1) = indsy(wav_propsy.offset_L+1); resend(2,res+1)   = indsy(end-wav_propsy.offset_R);
            resstart(3,res+1) = indsz(wav_propsz.offset_L+1); resend(3,res+1)   = indsz(end-wav_propsz.offset_R);
            
            sig_out(psiinds_x, indsy, indsz, :) = sig_in((endx-length(psiinds_x)+1):endx, 1:endy, 1:endz, :);
            sig_out(phiinds_x, psiinds_y, indsz, :) = sig_in(1:(endx-length(psiinds_x)), (endy-length(psiinds_y)+1):endy, 1:endz, :);
            sig_out(phiinds_x, phiinds_y, psiinds_z, :) = sig_in(1:(endx-length(psiinds_x)), 1:(endy-length(psiinds_y)), (endz-length(psiinds_z)+1):endx, :);
            
            endx = endx - length(psiinds_x); endy = endy - length(psiinds_y); endz = endz - length(psiinds_z);
            indsx = indsx((wav_propsx.offset_L+1):2:(end-wav_propsx.offset_R)); 
            indsy = indsy((wav_propsy.offset_L+1):2:(end-wav_propsy.offset_R));
            indsz = indsz((wav_propsz.offset_L+1):2:(end-wav_propsz.offset_R));
        end
        sig_out(indsx, indsy, indsz, :) = sig_in(1:endx, 1:endy, 1:endz, :);
    end
end
