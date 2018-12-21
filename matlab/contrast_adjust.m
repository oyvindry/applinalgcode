function Z=contrast_adjust(X, epsilon)
    Z = X/255; % Maps the pixel values to [0,1]
    Z = (log(Z+epsilon) - log(epsilon))/...
             (log(1+epsilon)-log(epsilon));
    Z = Z*255; % Maps the values back to [0,255]
end