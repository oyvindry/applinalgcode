test_orthogonality('db2')
test_orthogonality('db4')
test_dwt_different_sizes('cdf97')
test_dwt_different_sizes('spline53')
test_dwt_different_sizes('pwl0')
test_dwt_different_sizes('pwl2')
test_dwt_different_sizes('haar')
test_dwt_different_sizes('spline4.4')
test_simple_dwt2()
test_bd_db_van(4)

m = log2(64/(2*2));
test_bd('spline53', m)
test_bd('pwl2', m)
m = log2(64/(2*4));
test_bd('cdf97', m)
test_spline44()
    
function test_dwt_different_sizes(wave_name)
    disp('Testing the DWT on different input sizes')
    m = 4;

    disp(sprintf('Testing 2D with one channel: %s', wave_name))
    img = rand(64);
    img2 = img;
    img2 = dwt_impl(img2, wave_name, m, 'symm', 'none', 2);
    img2 = idwt_impl(img2, wave_name, m, 'symm', 'none', 2);
    diff = max(max(abs(img2-img)));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 2D with three channels: %s', wave_name))
    img = rand(64, 64, 3);
    img2 = img;
    img2 = dwt_impl(img2, wave_name, m);
    img2 = idwt_impl(img2, wave_name, m);
    diff = max(max(max(abs(img2-img))));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 1D with one channel: %s', wave_name))
    sd = rand(64,1);
    sd2 = sd;
    sd2 = dwt_impl(sd2, wave_name, m);
    sd2 = idwt_impl(sd2, wave_name, m);
    diff = max(abs(sd2-sd));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 1D with two channels: %s', wave_name))
    sd = rand(64,2);
    sd2 = sd;
    sd2 = dwt_impl(sd2, wave_name, m);
    sd2 = idwt_impl(sd2, wave_name, m);
    diff = max(max(abs(sd2-sd)));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 3D with one channel: %s', wave_name))
    sd = rand(64,64,64);
    sd2 = sd;
    sd2 = dwt_impl(sd2, wave_name, m, 'symm', 'none', 3);
    sd2 = idwt_impl(sd2, wave_name, m, 'symm', 'none', 3);
    diff = max(max(max(abs(sd2-sd))));
    assert(diff ~= 0 && diff < 1E-13)
    
    disp(sprintf('Testing 3D with three channels: %s', wave_name))
    sd = rand(64,64,64,3);
    sd2 = sd;
    sd2 = dwt_impl(sd2, wave_name, m);
    sd2 = idwt_impl(sd2, wave_name, m);
    diff = max(max(max(max(abs(sd2-sd)))));
    assert(diff ~= 0 && diff < 1E-13)
end
    
function test_orthogonality(wave_name)
    disp('Testing orthonormal wavelets:')
    x0 = rand(32,1);
    
    disp(sprintf('Testing that the IDWT inverts the DWT: %s', wave_name))
    x = x0;
    x = dwt_impl(x, wave_name, 2, 'per');
    x = idwt_impl(x, wave_name, 2, 'per');
    diff = max(abs(x-x0));
    assert(diff ~= 0 && diff < 1E-13)

    disp(sprintf('See that the wavelet transform equals the dual wavelet transform: %s', wave_name))
    x = x0;
    x = dwt_impl(x, wave_name, 2, 'per', 'none', 0, 1);
    x0 = dwt_impl(x0, wave_name, 2, 'per', 'none', 0, 0);
    diff = max(abs(x-x0));
    assert(diff ~= 0 && diff < 1E-13)

    disp(sprintf('Apply the transpose, to see that the transpose equals the inverse: %s', wave_name))
    x = x0;
    x = dwt_impl(x, wave_name, 2, 'per', 'none', 0, 0, 1);
    x = dwt_impl(x, wave_name, 2, 'per', 'none', 0, 0, 0);
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
    maxval = max(abs(x((64/2^m+2):65)));
    assert( maxval < 1E-11);
    x=idwt_impl(x,  wave_name, m, 'bd', 'bd_pre');
    diff = max(abs(res-x));
    assert(diff < 1E-12)
end
        
function test_bd_db_van(N)
    for k=1:N
        disp(sprintf('Testing bd db%i',k))
        m = floor(log2(64/(2*k+1))); % Max number of levels
        for s=0:(k-1)
            res = ((1:64).^s)';
            x = dwt_impl (res, sprintf('db%i',k),  m, 'bd', 'bd_pre');
            maxval = max(abs(x((64/2^m+1):64)));
            assert( maxval < 1E-8)
            x = idwt_impl(x,   sprintf('db%i',k), m,  'bd', 'bd_pre');
            diff = max(abs(res-x));
            assert(diff ~= 0 && diff < 1E-8)
        end
    end
end

function test_spline44()
    disp('Testing spline4.4')
    m = 2;
    for s=0:3
        res = ((1:65).^s)';
        x=dwt_impl(res, 'spline4.4',  m, 'bd', 'bd_pre');
        maxval = max(abs(x((64/2^m+2):65)));
        assert( maxval < 1E-7)
        x=idwt_impl(x, 'spline4.4',  m, 'bd', 'bd_pre');
        diff = max(abs(res-x));
        assert(diff ~= 0 && diff < 1E-7)
    end
end

% x = double(imread('images/lena.png'));
% x=dwt_impl(x, 'cdf97', 1, 'bd', 'bd_pre'); 
% x=idwt_impl(x, 'cdf97', 1, 'bd', 'bd_pre');
% imshow(uint8(x))