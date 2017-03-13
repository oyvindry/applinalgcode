function [img1,img2,img3]=mmsubbands(m)
    img = CreateExcerpt();
    [l1, l2, l3] = size(img);
    X = zeros(3,l1,l2,l3);
    
    img1 = img;
    img1 = DWT2Impl(img1,m,'cdf53');
    img1((l1/2^(m-1)+1):end, :, :) = 0;
    img1(1:(l1/2^(m-1)), (l2/2^(m-1)+1):end, :) = 0;
    img1((l1/2^m+1):(l1/2^(m-1)), (l2/2^m+1):(l2/2^(m-1)), :) = 0;
    img2 = img1;
    img1 = IDWT2Impl(img1, m, 'cdf53');
    
    img2((l1/2^m+1):(l1/2^(m-1)), 1:(l2/2^m), :) = 0;
    img3 = img2;
    img2 = IDWT2Impl(img2, m, 'cdf53');
    
    img3(1:(l1/2^m), (l2/2^m+1):(l2/2^(m-1)), :) = 0;
    img3 = IDWT2Impl(img3, m, 'cdf53');
    
    
    img1 = mapto01(img1);
    img1 = img1*255;
    img2 = mapto01(img2);
    img2 = img2*255;
    img3 = mapto01(img3);
    img3 = img3*255;
end