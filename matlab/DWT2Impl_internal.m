function x = DWT2Impl_internal(x, m, fx, fy, bd_mode, prefilterx, prefiltery, wav_propsx, wav_propsy, data_layout)
    lastdim = 1;
    if length(size(x)) == 3
        lastdim = size(x, 3);
    end 
    indsx = 1:size(x,1); indsy = 1:size(x,2);

    % preconditioning   
    x=tensor2_kernel(x, indsx, indsy, @(x,bd_mode) prefilterx(x, 1), @(x,bd_mode) prefiltery(x, 1), lastdim, bd_mode);
         
    for res = 0:(m - 1)
        x=tensor2_kernel(x, indsx, indsy, fx, fy, lastdim, bd_mode);  
        indsx = indsx((wav_propsx.offset_L+1):2:(end-wav_propsx.offset_R));
        indsy = indsy((wav_propsy.offset_L+1):2:(end-wav_propsy.offset_R));
    end
    
    % postconditioning
    x=tensor2_kernel(x, indsx, indsy, @(x,bd_mode) prefilterx(x, 0), @(x,bd_mode) prefiltery(x, 0), lastdim, bd_mode);  
    
    x = reorganize_coeffs2_forward(x, m, wav_propsx, wav_propsy, data_layout);   
end     
        
function sig_out=reorganize_coeffs2_forward(sig_in, m, wav_propsx, wav_propsy, data_layout)
    sig_out = sig_in;
    if strcmpi(data_layout, 'resolution')
        indsx = 1:size(sig_in,1); indsy = 1:size(sig_in,2);
        endx = size(sig_in,1); endy = size(sig_in,2); 
        for res=1:m
            psiinds_x = [indsx(1:wav_propsx.offset_L) indsx((wav_propsx.offset_L + 2):2:(end-wav_propsx.offset_R)) indsx((end-wav_propsx.offset_R+1):end)]; % psi-indices
            psiinds_y = [indsy(1:wav_propsy.offset_L) indsy((wav_propsy.offset_L + 2):2:(end-wav_propsy.offset_R)) indsy((end-wav_propsy.offset_R+1):end)];
            phiinds_x = indsx((wav_propsx.offset_L + 1):2:(end-wav_propsx.offset_R));
            
            sig_out( (endx-length(psiinds_x)+1):endx, 1:endy, :) = sig_in(psiinds_x,indsy,:);
            sig_out( 1:(endx-length(psiinds_x)), (endy-length(psiinds_y)+1):endy, :) = sig_in(phiinds_x,psiinds_y,:);
            
            endx = endx - length(psiinds_x); endy = endy - length(psiinds_y);
            indsx = indsx((wav_propsx.offset_L+1):2:(end-wav_propsx.offset_R)); 
            indsy = indsy((wav_propsy.offset_L+1):2:(end-wav_propsy.offset_R));
        end
        sig_out(1:endx, 1:endy, :) = sig_in(indsx, indsy, :);
    end
end