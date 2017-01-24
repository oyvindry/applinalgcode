function x = IDWTKernelOrtho(x, filters, mode, dual)
    % For usage see - DWTKernelOrtho
    
    N = size(x, 1);
    y1 = 0; y2 = 0;
    
    if mode >= 2
       y1 = filters.AL*x(1:size(filters.AL,2));
       y2 = filters.AR*x((N-size(filters.AR,2)+1):N);
    end
    if dual
        stepnr = 1;
        if mod(size(filters.lambdas, 1), 2) == 1
            x = liftingstepodd(-filters.lambdas(stepnr, 2), -filters.lambdas(stepnr, 1), x, mode);
            stepnr = stepnr + 1;
        end
  
        while stepnr < size(filters.lambdas, 1)
            x = liftingstepeven(-filters.lambdas(stepnr, 2), -filters.lambdas(stepnr, 1), x, mode);
            stepnr = stepnr + 1;
            x = liftingstepodd(-filters.lambdas(stepnr, 2), -filters.lambdas(stepnr, 1), x, mode);
            stepnr = stepnr + 1;
        end
 
        x(1:2:N, :) = x(1:2:N, :)*filters.alpha;
        x(2:2:N, :) = x(2:2:N, :)*filters.beta;
    else
        stepnr = 1;
        if mod(size(filters.lambdas, 1), 2) == 1
            x = liftingstepeven(filters.lambdas(stepnr, 1), filters.lambdas(stepnr, 2), x, mode);
            stepnr = stepnr + 1;
        end
  
        while stepnr < size(filters.lambdas, 1)
            x = liftingstepodd(filters.lambdas(stepnr, 1), filters.lambdas(stepnr, 2), x, mode);
            stepnr = stepnr + 1;
            x = liftingstepeven(filters.lambdas(stepnr, 1), filters.lambdas(stepnr, 2), x, mode);
            stepnr = stepnr + 1;
        end
 
        x(1:2:N, :)=x(1:2:N, :)/filters.alpha;
        x(2:2:N, :)=x(2:2:N, :)/filters.beta;
    end
    if mode >= 2
        x(1:size(filters.AL,1)) = x(1:size(filters.AL,1)) + y1;
        x((N-size(filters.AR,1)+1):N) = x((N-size(filters.AR,1)+1):N) + y2;
        if mode == 3
            x(1:size(filters.A_L_pre_inv,1)) = filters.A_L_pre_inv*x(1:size(filters.A_L_pre_inv,1));
            x((end-size(filters.A_R_pre_inv,1)+1):end) = filters.A_R_pre_inv*x((end-size(filters.A_R_pre_inv,1)+1):end);
        end
    end
end
