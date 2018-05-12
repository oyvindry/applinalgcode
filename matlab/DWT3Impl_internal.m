function x = DWT3Impl_internal(x, m, fx, fy, fz, bd_mode, prefilterx, prefiltery, prefilterz, wav_propsx, wav_propsy, wav_propsz, data_layout)
    lastdim = 1;
    if length(size(x)) == 4
        lastdim = size(x, 4);
    end 
    indsx = 1:size(x,1); indsy = 1:size(x,2); indsz = 1:size(x,3);

    % preconditioning   
    x = tensor3_kernel(x, indsx, indsy, indsz, @(x,bd_mode) prefilterx(x, 1), @(x,bd_mode) prefiltery(x, 1), @(x,bd_mode) prefilterz(x, 1), lastdim, bd_mode);
         
    for res = 0:(m - 1)
        x=tensor3_kernel(x, indsx, indsy, indsz, fx, fy, fz, lastdim, bd_mode);  
        indsx = indsx((wav_propsx.offset_L+1):2:(end-wav_propsx.offset_R));
        indsy = indsy((wav_propsy.offset_L+1):2:(end-wav_propsy.offset_R));
        indsz = indsz((wav_propsz.offset_L+1):2:(end-wav_propsz.offset_R));
    end
    
    % postconditioning
    x=tensor3_kernel(x, indsx, indsy, indsz, @(x,bd_mode) prefilterx(x, 0), @(x,bd_mode) prefiltery(x, 0), @(x,bd_mode) prefilterz(x, 0), lastdim, bd_mode);  
    
    x = reorganize_coeffs3_forward(x, m, wav_propsx, wav_propsy, wav_propsz, data_layout);   
end

function sig_out=reorganize_coeffs3_forward(sig_in, m, wav_propsx, wav_propsy, wav_propsz, data_layout)
    sig_out = sig_in;
    if strcmpi(data_layout, 'resolution')
        indsx = 1:size(sig_in,1); indsy = 1:size(sig_in,2); indsz = 1:size(sig_in,3);
        endx = size(sig_in,1); endy = size(sig_in,2); endz = size(sig_in,3); 
        for res=1:m
            psiinds_x = [indsx(1:wav_propsx.offset_L) indsx((wav_propsx.offset_L + 2):2:(end-wav_propsx.offset_R)) indsx((end-wav_propsx.offset_R+1):end)]; % psi-indices
            psiinds_y = [indsy(1:wav_propsy.offset_L) indsy((wav_propsy.offset_L + 2):2:(end-wav_propsy.offset_R)) indsy((end-wav_propsy.offset_R+1):end)];
            psiinds_z = [indsz(1:wav_propsz.offset_L) indsz((wav_propsz.offset_L + 2):2:(end-wav_propsz.offset_R)) indsz((end-wav_propsz.offset_R+1):end)];
            phiinds_x = indsx((wav_propsx.offset_L + 1):2:(end-wav_propsx.offset_R));
            phiinds_y = indsy((wav_propsy.offset_L + 1):2:(end-wav_propsy.offset_R));
            
            sig_out( (endx-length(psiinds_x)+1):endx, 1:endy, 1:endz, :) = sig_in(psiinds_x, indsy, indsz, :);
            sig_out( 1:(endx-length(psiinds_x)), (endy-length(psiinds_y)+1):endy, 1:endz, :) = sig_in(phiinds_x, psiinds_y, indsz, :);
            sig_out( 1:(endx-length(psiinds_x)), 1:(endy-length(psiinds_y)), (endz-length(psiinds_z)+1):endx, :) = sig_in(phiinds_x, phiinds_y, psiinds_z, :);
            
            endx = endx - length(psiinds_x); endy = endy - length(psiinds_y); endz = endz - length(psiinds_z);
            indsx = indsx((wav_propsx.offset_L+1):2:(end-wav_propsx.offset_R)); 
            indsy = indsy((wav_propsy.offset_L+1):2:(end-wav_propsy.offset_R));
            indsz = indsz((wav_propsz.offset_L+1):2:(end-wav_propsz.offset_R));
        end
        sig_out(1:endx, 1:endy, 1:endz, :) = sig_in(indsx, indsy, indsz, :);
    end
end