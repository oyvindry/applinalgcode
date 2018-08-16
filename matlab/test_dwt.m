test_dwt_different_sizes()
test_kernel_ortho()
test_orthogonality()
test_kernel('cdf97')
test_kernel('cdf53')
test_kernel('pwl0')
test_kernel('pwl2')
test_kernel('haar')
test_kernel('spline4.4')
test_simple_dwt2()
test_bd_db_van(4)
test_bd('cdf53', 4)
test_bd('pwl2', 4)
test_spline44()
    
function test_kernel(wave_name)
    disp(sprintf('Testing %s, 1D', wave_name))
    % res = rand(16,1);
    res = [1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16];
    x = res;
    x = dwt_impl(x, wave_name, 2);
    x = idwt_impl(x, wave_name, 2);
    diff = max(abs(x-res));
    assert(diff < 1E-13)
    
    disp(sprintf('Testing %s, 1D, two channels', wave_name))
    res = rand(16,2);
    x = res;
    x = dwt_impl(x, wave_name, 2);
    x = idwt_impl(x, wave_name, 2);
    diff = max(max(abs(x-res)));
    assert(diff < 1E-13)
    
    disp(sprintf('Testing %s, 2D', wave_name))
    res = rand(16,16);
    x = res;
    x = dwt_impl(x, wave_name, 2, 'symm', 'none', 2);
    x = idwt_impl(x, wave_name, 2, 'symm', 'none', 2);
    diff = max(max(max(abs(x-res))));
    assert(diff < 1E-13)
    
    disp(sprintf('Testing %s, 2D, 3 channels', wave_name))
    res = rand(16,16,3);
    x = res;
    x = dwt_impl(x, wave_name, 2);
    x = idwt_impl(x, wave_name, 2);
    diff = max(max(max(abs(x-res))));
    assert(diff < 1E-13)
end
    
function test_kernel_ortho()
    disp('Testing orthonormal wavelets')
    res = rand(16,1); % only this assumes that N is even
    
    disp('Testing that the reverse inverts the forward transform')
    x = res;
    x = dwt_impl(x, 'db4', 2);
    x = idwt_impl(x, 'db4', 2);
    diff = max(abs(x-res));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp('Testing that the transform is orthogonal, i.e. that the transform and its dual are equal')
    x = res;
    x = dwt_impl(x, 'db4', 2, 'per');
    res = dwt_impl(res, 'db4', 2, 'per', 'none', 0, 1);
    diff = max(abs(x-res));
    assert(diff ~= 0 && diff < 1E-13)
end
    
function test_dwt_different_sizes()
    disp('Testing the DWT on different input sizes')
    m = 4;

    disp('Testing the DWT for greyscale image')
    img = rand(32);
    img2 = img;
    img2 = dwt_impl(img2, 'cdf97', m, 'symm', 'none', 2);
    img2 = idwt_impl(img2, 'cdf97', m, 'symm', 'none', 2);
    diff = max(max(abs(img2-img)));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp('Testing the DWT for RGB image')
    img = rand(32, 32, 3);
    img2 = img;
    img2 = dwt_impl(img2, 'cdf97', m);
    img2 = idwt_impl(img2, 'cdf97', m);
    diff = max(max(max(abs(img2-img))));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp('Testing the DWT for sound with one channel')
    sd = rand(32,1);
    sd2 = sd;
    sd2 = dwt_impl(sd2, 'cdf97', m);
    sd2 = idwt_impl(sd2, 'cdf97', m);
    diff = max(abs(sd2-sd));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp('Testing the DWT for sound with two channels')
    sd = rand(32,2);
    sd2 = sd;
    sd2 = dwt_impl(sd2, 'cdf97', m);
    sd2 = idwt_impl(sd2, 'cdf97', m);
    diff = max(max(abs(sd2-sd)));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp('Testing 3D with one channel')
    sd = rand(32,32,32);
    sd2 = sd;
    sd2 = dwt_impl(sd2, 'cdf97', m, 'symm', 'none', 3);
    sd2 = idwt_impl(sd2, 'cdf97', m, 'symm', 'none', 3);
    diff = max(max(max(abs(sd2-sd))));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp('Testing 3D with two channels')
    sd = rand(32,32,32,3);
    sd2 = sd;
    sd2 = dwt_impl(sd2, 'cdf97', m);
    sd2 = idwt_impl(sd2, 'cdf97', m);
    diff = max(max(max(max(abs(sd2-sd)))));
    assert(diff ~= 0 && diff < 1E-13)
end
    
function test_orthogonality()
    disp('Testing that the wavelet and the dual wavelet are equal for orthonormal wavelets')
    x0 = rand(32,1);
    
    disp('Testing that the IDWT inverts the DWT')
    x = x0;
    x = dwt_impl(x, 'db4', 2, 'per');
    x = idwt_impl(x, 'db4', 2, 'per');
    diff = max(abs(x-x0));
    assert(diff ~= 0 && diff < 1E-13)

    disp('See that the wavelet transform equals the dual wavelet transform')
    x = x0;
    x = dwt_impl(x, 'db4', 2, 'per', 'none', 0, 1);
    x0 = dwt_impl(x0, 'db4', 2, 'per', 'none', 0, 0);
    diff = max(abs(x-x0));
    assert(diff ~= 0 && diff < 1E-13)

    disp('Apply the transpose, to see that the transpose equals the inverse')
    x = x0;
    x = dwt_impl(x, 'db4', 2, 'per', 'none', 0, 0, 1);
    x = dwt_impl(x, 'db4', 2, 'per', 'none', 0, 0, 0);
    diff = max(abs(x-x0));
    assert(diff ~= 0 && diff < 1E-13)
end

function test_simple_dwt2()
    disp('Testing simple DWT2')
    img = rand(32, 32, 3);
    img2 = img;
% Begin simple_dwt2
    f = @(x, bd_mode) dwt_impl(x, 'cdf97', 4, bd_mode, 'none', 1);
    img = tensor2_impl(img, f, f, 'symm');
% End simple_dwt2
% Begin simple_idwt2
    invf = @(x, bd_mode) idwt_impl(x, 'cdf97', 4, bd_mode, 'none', 1);
    img = tensor2_impl(img, invf, invf, 'symm');
% End simple_idwt2
    diff = max(max(max(abs(img2-img))));
    assert(diff ~= 0 && diff < 1E-13)
end

function test_bd(wave_name, m)
    disp(sprintf('Testing bd %s', wave_name))
    res = (1:65)';
    x=dwt_impl(res, wave_name,  m, 'bd', 'bd_pre');
    maxval = max(abs(x(6:64)));
    assert( maxval < 1E-13)
    x=idwt_impl(x,  wave_name, m, 'bd', 'bd_pre');
    diff = max(abs(res-x));
    assert(diff < 1E-13)
end
        
function test_bd_db_van(N)
    for k=1:N
        disp(sprintf('Testing bd db%i',k))
        for s=0:(k-1)
            res = ((1:64).^s)';
            x = dwt_impl (res, sprintf('db%i',k),  2, 'bd', 'bd_pre');
            maxval = max(abs(x(17:64)));
            assert( maxval < 1E-8)
            x = idwt_impl(x,   sprintf('db%i',k), 2,  'bd', 'bd_pre');
            diff = max(abs(res-x));
            assert(diff ~= 0 && diff < 1E-8)
        end
    end
end

function test_spline44()
    disp('Testing spline4.4')
    for s=0:3
        res = (1:63)';
        x=dwt_impl(res, 'spline4.4',  1, 'bd', 'bd_pre');
        maxval = max(abs(x(32:63)));
        assert( maxval < 1E-8)
        x=idwt_impl(x, 'spline4.4',  1, 'bd', 'bd_pre');
        diff = max(abs(res-x));
        assert(diff ~= 0 && diff < 1E-8)
    end
end

% x = double(imread('images/lena.png'));
% x=dwt_impl(x, 'cdf97', 1, 'bd', 'bd_pre'); 
% x=idwt_impl(x, 'cdf97', 1, 'bd', 'bd_pre');
% imshow(uint8(x))