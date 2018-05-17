function [img1,img2,img3]=mmsubbands(m)
    img = create_excerpt();
    [l1, l2, l3] = size(img);
    X = zeros(3,l1,l2,l3);
    
    img1 = img;
    img1 = dwt_impl(img1, 'cdf53', m);
    img1((l1/2^(m-1)+1):end, :, :) = 0;
    img1(1:(l1/2^(m-1)), (l2/2^(m-1)+1):end, :) = 0;
    img1((l1/2^m+1):(l1/2^(m-1)), (l2/2^m+1):(l2/2^(m-1)), :) = 0;
    img2 = img1;
    img1 = idwt_impl(img1, 'cdf53', m);
    
    img2((l1/2^m+1):(l1/2^(m-1)), 1:(l2/2^m), :) = 0;
    img3 = img2;
    img2 = idwt_impl(img2, 'cdf53', m);
    
    img3(1:(l1/2^m), (l2/2^m+1):(l2/2^(m-1)), :) = 0;
    img3 = idwt_impl(img3, 'cdf53', m);
    
    
    img1 = mapto01(img1);
    img1 = img1*255;
    img2 = mapto01(img2);
    img2 = img2*255;
    img3 = mapto01(img3);
    img3 = img3*255;
end