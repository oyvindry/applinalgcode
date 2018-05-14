function x=tensor2_kernel(x, indsx, indsy, fx, fy, lastdim, bd_mode)
    for col = indsy
        Y2 = zeros(length(indsx), lastdim);
        Y2(:, :) = x(indsx, col, :);
        x(indsx, col, :) = fx(Y2, bd_mode);
    end
        
    for row = indsx
        Y1=zeros(length(indsy), lastdim); 
        Y1(:, :) = x(row, indsy, :);
        x(row, indsy, :) = fy(Y1, bd_mode);
    end
end