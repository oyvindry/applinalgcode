function cascade_alg(m, a, b, wave_name, scaling, dual)
    coords = zeros((b-a)*2^m, 1);
    if scaling
        coords(1) = 1;
    else
        coords(b - a + 1) = 1;
    end
    t = linspace(a, b, (b-a)*2^m);
    coords = dwt_impl(coords, m, wave_name, 0, 1, 0, dual);
    coords = [ coords((b*2^m+1):((b-a)*2^m)); coords(1:(b*2^m)) ];
    figure()
    plot(t, 2^(m/2)*coords, 'k-')