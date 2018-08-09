function data=tensor3_impl(data, fx, fy, fz, lastdim, bd_mode)
    [xlen, ylen, zlen]=size(data);
    for y = 1:ylen
        for z = 1:zlen
            data_temp = zeros(xlen, lastdim);
            data_temp(:, :) = data(:, y, z, :);
            data(:, y, z, :) = fx(data_temp, bd_mode);
        end
    end
        
    for z = 1:zlen
        for x = 1:xlen
            data_temp = zeros(ylen, lastdim);
            data_temp(:, :) = data(x, :, z, :);
            data(x, :, z, :) = fy(data_temp, bd_mode);
        end
    end
    
    for x = 1:xlen
        for y = 1:ylen
            data_temp = zeros(zlen, lastdim);
            data_temp(:, :) = data(x, y, :, :);
            data(x, y, :, :) = fz(data_temp, bd_mode);
        end
    end
end