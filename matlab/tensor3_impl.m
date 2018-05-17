function data=tensor3_impl(data, indsx, indsy, indz, fx, fy, fz, lastdim, bd_mode)
    for y = indsy
        for z = indsz
            data_temp = zeros(length(indsx), lastdim);
            data_temp(:, :) = data(indsx, y, z, :);
            data(indsx, y, z, :) = fx(data_temp, bd_mode);
        end
    end
        
    for z = indsz
        for x = indsx
            data_temp = zeros(length(indsy), lastdim);
            data_temp(:, :) = data(x, indsy, z, :);
            data(x, indsy, z, :) = fy(data_temp, bd_mode);
        end
    end
    
    for x = indsx
        for y = indsy
            data_temp = zeros(length(indsz), lastdim);
            data_temp(:, :) = data(x, y, indsz, :);
            data(x, y, indsz, :) = fz(data_temp, bd_mode);
        end
    end
end