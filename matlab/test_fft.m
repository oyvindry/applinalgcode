testfft()
testdct()
    
function testfft()
    x1 = rand(32,1);
    x2 = rand(32,3);
    x2copy = x2;
    
    disp('Testing dft_impl')
    x = fft(x1);
    x0 = dft_impl(x1);
    diff = max(abs(x0-x));
    assert(diff < 1E-13)
    
    x = fft(x2);
    x0 = dft_impl(x2);
    diff = max(max(abs(x0-x)));
    assert(diff < 1E-13)
    
    disp('Testing fft_kernel_standard')
    x = fft(x1);
    x1copy = x1;
    x1copy = fft_impl(x1copy, @fft_kernel_standard);
    diff = max(abs(x1copy-x));
    assert(diff < 1E-13)
    
    x = fft(x2);
    x2copy = x2;
    x2copy = fft_impl(x2copy, @fft_kernel_standard);
    diff = max(max(abs(x2copy-x)));
    assert(diff < 1E-13)
    
    disp('Testing fft_kernel_nonrec')
    x = fft(x1);
    x1copy = x1;
    x1copy = fft_impl(x1copy, @fft_kernel_nonrec);
    diff = max(abs(x1copy-x));
    assert(diff < 1E-13)
    
    x = fft(x2);
    x2copy = x2;
    x2copy = fft_impl(x2copy, @fft_kernel_nonrec);
    diff = max(max(abs(x2copy-x)));
    assert(diff < 1E-13)
    
    disp('Testing fft_kernel_splitradix')
    x = fft(x1);
    x1copy = x1;
    x1copy = fft_impl(x1copy, @fft_kernel_splitradix);
    diff = max(abs(x1copy-x));
    assert(diff < 1E-13)
    
    x = fft(x2);
    x2copy = x2;
    x2copy = fft_impl(x2copy, @fft_kernel_splitradix);
    diff = max(max(abs(x2copy-x)));
    assert(diff < 1E-13)
    
    disp('Testing that fft_kernel_standard(reverse) inverts fft_kernel_standard(forward)') 
    x1copy = x1;
    x1copy = fft_impl(x1copy, @fft_kernel_standard);
    x1copy = fft_impl(x1copy, @fft_kernel_standard, 0);
    diff = max(abs(x1copy-x1));
    assert(diff < 1E-13)
    
    x2copy = x2;
    x2copy = fft_impl(x2copy, @fft_kernel_standard);
    x2copy = fft_impl(x2copy, @fft_kernel_standard, 0); 
    diff = max(max(abs(x2copy-x2)));
    assert(diff < 1E-13)
    
    disp('Testing that fft_kernel_nonrec(reverse) inverts fft_kernel_nonrec(forward)') 
    x1copy = x1;
    x1copy = fft_impl(x1copy, @fft_kernel_nonrec);
    x1copy = fft_impl(x1copy, @fft_kernel_nonrec, 0); 
    diff = max(abs(x1copy-x1));
    assert(diff < 1E-13)
    
    x2copy = x2;
    x2copy = fft_impl(x2copy, @fft_kernel_nonrec);
    x2copy = fft_impl(x2copy, @fft_kernel_nonrec, 0); 
    diff = max(max(abs(x2copy-x2)));
    assert(diff < 1E-13)
    
    disp('Testing that fft_kernel_splitradix(reverse) inverts fft_kernel_splitradix(forward)') 
    x1copy = x1;
    x1copy = fft_impl(x1copy, @fft_kernel_splitradix);
    x1copy = fft_impl(x1copy, @fft_kernel_splitradix, 0); 
    diff = max(abs(x1copy-x1));
    assert(diff < 1E-13)
    
    x2copy = x2;
    x2copy = fft_impl(x2copy, @fft_kernel_splitradix);
    x2copy = fft_impl(x2copy, @fft_kernel_splitradix, 0); 
    diff = max(max(abs(x2copy-x2)));
    assert(diff < 1E-13)
    
    disp('Testing that dft_impl(reverse) inverts dft_impl(forward)')
    x0 = dft_impl(dft_impl(x1), 0);
    diff = max(abs(x0-x1));
    assert(diff < 1E-13)
    
    x0 = dft_impl(dft_impl(x2), 0);
    diff = max(max(abs(x0-x2)));
    assert(diff < 1E-13)
end

function testdct()
    x1 = rand(32,1);
    x2 = rand(32,3);
    x2copy = x2;
    
    disp('Testing dct_impl')
    x = dct(x1);
    x1copy = x1;
    x1copy = dct_impl(x1copy);
    diff = max(abs(x1copy-x));
    assert(diff < 1E-13)
    
    x = dct(x2);
    x2copy = x2;
    x2copy = dct_impl(x2copy);
    diff = max(max(abs(x2copy-x)));
    assert(diff < 1E-13)
    
    
    disp('Testing that idct_impl inverts dct_impl')
    x1copy = x1;
    x1copy = dct_impl(x1copy);
    x1copy = idct_impl(x1copy); 
    diff = max(abs(x1copy-x1));
    assert(diff < 1E-13)
    
    x2copy = x2;
    x2copy = dct_impl(x2copy);
    x2copy = idct_impl(x2copy); 
    diff = max(max(abs(x2copy-x2)));
    assert(diff < 1E-13)
end