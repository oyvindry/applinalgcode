function J0 = computeJ0(vm)
    J0 = 0;
    while (2^J0 <= 2*vm) % This is not optimal in all cases, but it is 
                         % necessary in order to make a unified interface with 
                         % the non-boundary wavelets. Do not change it! 
        J0 = J0 + 1;
    end
end

