function [wav_props, dual_wav_props=find_wavprops(m, wave_name, bd_mode, lengthsignal)
    wav_props.wave_name = wave_name; dual_wav_props.wave_name = wave_name;
    wav_props.m = m; dual_wav_props.m = m;
    wav_props.lengthsignal = lengthsignal; dual_wav_props.lengthsignal = lengthsignal;
    
    if (strcmpi(wave_name(1:2), 'db'))
        N = str2double(wave_name(3:end));
        if N > 1 
            [wav_props,dual_wav_props] = wav_props_ortho(N, wav_props,dual_wav_props, bd_mode);
        end
    elseif (strcmpi(wave_name(1:3), 'sym'))
        N = str2double(wave_name(4:end));
        if N > 1
            [wav_props,dual_wav_props] = wav_props_ortho(N, wav_props,dual_wav_props, bd_mode, 1);
        end
    elseif strcmpi(wave_name(1:6), 'spline')
        N = str2double(wave_name(7));
        Ntilde = str2double(wave_name(9));
        [wav_props, dual_wav_props] = wav_props_biortho(N, Ntilde, wav_props, dual_wav_props, bd_mode);
    end
end

function [wav_props, dual_wav_props]=wav_props_ortho(N, wav_props, dual_wav_props, bd_mode, type)
    % N:    Number of vanishing moments
    % type: The type of orthonormal wavelet.
    %       0: Daubechies wavelets with minimum phase (default)  
    %       1: Symmlets - wavelets with close to linear phase (almost symmetric)
    % 
    if (nargin == 3)
        type = 0;
    end
        
    if (type == 0)
        [wav_props.h0, wav_props.h1, wav_props.g0, wav_props.g1] = h0h1computeortho(N);
    elseif (type == 1)
        [wav_props.h0, wav_props.h1, wav_props.g0, wav_props.g1] = h0h1computesym(N);
    end
    dual_wav_props.h0 = flip(wav_props.g0); dual_wav_props.h1 = flip(wav_props.g1); dual_wav_props.g0 = flip(wav_props.h0); dual_wav_props.g1 = flip(wav_props.h1);
    [wav_props.lambdas, wav_props.alpha, wav_props.beta, wav_props.last_even] = liftingstepscomputeortho(wav_props.h0, wav_props.h1);
    dual_wav_props.lambdas = -fliplr(wav_props.lambdas); dual_wav_props.alpha = 1/wav_props.alpha; dual_wav_props.beta = 1/wav_props.beta;
    dual_wav_props.last_even = ~wav_props.last_even;
    
    if strcmpi(bd_mode, 'bd')
        [wav_props, dual_wav_props] = wav_props_ortho_bd(N, wav_props, dual_wav_props);
    end 
end
    

function [wav_props, dual_wav_props]=wav_props_ortho_bd(N, wav_props, dual_wav_props)
    % Compute K_L and K_R
    K_L = N; K_R = N;
    s = mod(wav_props.lengthsignal, 2^wav_props.m);
    if s > 0
        toadd = 2^wav_props.m - s;
        K_L = K_L + floor(toadd/2);
        K_R = K_R + ceil(toadd/2);
    end
    wav_props.offset_L = K_L - N; wav_props.offset_R = K_R - N; dual_wav_props.offset_L = K_L - N; dual_wav_props.offset_R = K_R - N; 
        
    % First the right edge
    wav_props.h0 = flip(wav_props.h0); wav_props.h1 = flip(wav_props.h1); wav_props.g0 = flip(wav_props.g0); wav_props.g1 = flip(wav_props.g1);
    [W, A_pre, A_pre_inv] = bw_compute_left_ortho(wav_props.g0, wav_props.g1, N, K_R); % Lower right (3N-1)x(2N) matrix
    wav_props.A_R_pre = fliplr(flipud(A_pre)); wav_props.A_R_pre_inv = fliplr(flipud(A_pre_inv));
    dual_wav_props.A_R_pre = wav_props.A_R_pre; dual_wav_props.A_R_pre_inv = wav_props.A_R_pre_inv;
        
    WR = zeros(size(W));
    for k = (size(W,2)-(K_R-N)):(-2):1
        WR(:,[k-1 k]) = W( size(W,1):(-1):1, size(W,2) + 1 - [k k-1]); 
    end
       
    % Then the left edge
    wav_props.h0 = flip(wav_props.h0); wav_props.h1 = flip(wav_props.h1); wav_props.g0 = flip(wav_props.g0); wav_props.g1 = flip(wav_props.g1);
    [WL, wav_props.A_L_pre, wav_props.A_L_pre_inv] = bw_compute_left_ortho(wav_props.g0, wav_props.g1, N, K_L); % Upper left (3N-1)x(2N) matrix
    dual_wav_props.A_L_pre = wav_props.A_L_pre; dual_wav_props.A_L_pre_inv = wav_props.A_L_pre_inv;
    
    % Compute the left and right parts of the IDWT
    M = 6*N;
    x = IDWTImpl(eye(M), 1, wav_props.wave_name, 'none', 'none', 0, 0, 'time');
        
    [w1, w2] = size(WL);
    wav_props.A_L = WL - x(1:w1,1:w2); dual_wav_props.A_L = wav_props.A_L;
    wav_props.A_R = WR - x((M-w1+1):M,(M-w2+1):M); dual_wav_props.A_R = wav_props.A_R;
end


function [wav_props, dual_wav_props] = wav_props_biortho(N, Ntilde, wav_props, dual_wav_props, bd_mode)
    [wav_props.h0, wav_props.h1, wav_props.g0, wav_props.g1]=compute_spline_filters(N, Ntilde);
    dual_wav_props.g0 = wav_props.h0;
    dual_wav_props.g1 = wav_props.h1;
    dual_wav_props.h0 = wav_props.g0; 
    dual_wav_props.h1 = wav_props.g1;
    if strcmpi(bd_mode, 'bd')
        [wav_props, dual_wav_props]=wav_props_biortho_bd(N, Ntilde, wav_props, dual_wav_props);
    end 
end   

function [wav_props, dual_wav_props]=wav_props_biortho_bd(N, Ntilde, wav_props, dual_wav_props)
    Nprime = max(N,Ntilde);
    R = (length(wav_props.g0)-1)/2; Rtilde = (length(wav_props.h0)-1)/2;
    L  = -R; Ltilde = -Rtilde;
        
    K_L = -L; K_L_tilde = K_L + Ntilde - N;
    if K_L_tilde < -Ltilde
        K_L_tilde = -Ltilde;
        K_L = K_L_tilde + N - Ntilde;
    end
    K_R = K_L; K_R_tilde = K_L_tilde;
        
    s = wav_props.lengthsignal + R + L - 2*Nprime + K_L + K_R - 1;
    if mod(s,2^wav_props.m) > 0
        toadd = 2^wav_props.m - mod(s,2^wav_props.m);
        K_L = K_L + floor(toadd/2); K_L_tilde = K_L_tilde + floor(toadd/2);
        K_R = K_R + ceil(toadd/2);  K_R_tilde = K_R_tilde + ceil(toadd/2);
    end
    wav_props.offset_L = K_L - N; wav_props.offset_R = K_R - N; dual_wav_props.offset_L = wav_props.offset_L; dual_wav_props.offset_R = wav_props.offset_R;
        
        
    [WL, WLtilde, wav_props.A_L_pre, wav_props.A_L_pre_inv, dual_wav_props.A_L_pre, dual_wav_props.A_L_pre_inv] = bw_compute_left_biortho(wav_props.g0, wav_props.g1, N, wav_props.h0, wav_props.h1, Ntilde, K_L, K_L_tilde);
    [WR, WRtilde, wav_props.A_R_pre, wav_props.A_R_pre_inv, dual_wav_props.A_R_pre, dual_wav_props.A_R_pre_inv] = bw_compute_left_biortho(wav_props.g0, wav_props.g1, N, wav_props.h0, wav_props.h1, Ntilde, K_R, K_R_tilde);
    
    % Mirror WR and WRtilde
    W = zeros(size(WR)); 
    for k = (size(WR,2)-(K-N)):(-2):1
        W(:,[k-1 k]) = WR( size(WR,1):(-1):1, size(WR,2) + 1 - [k k-1]); 
    end
    WR = W;
    WRtilde = zeros(size(WRtilde));
    for k = (size(WRtilde,2)-(K-N)):(-2):1
        W(:,[k-1 k]) = WRtilde( size(WRtilde,1):(-1):1, size(WRtilde,2) + 1 - [k k-1]); 
    end
    WRtilde = W;
    
    M = max( [size(WL)+size(WR) size(WLtilde)+size(WRtilde)] );
    
    x1 = IDWTImpl(eye(M), 1, wav_props.wave_name, 'none', 'none', 0, 0, 'time');
    [w1, w2] = size(WL); wav_props.A_L = WL - x1(1:w1,1:w2);
    [w1, w2] = size(WR); wav_props.A_R = WR - x1((M-w1+1):M,(M-w2+1):M);
    
    x2 = IDWTImpl(eye(M), 1, wav_props.wave_name, 'none', 'none', 1, 0, 'time');
    [w1, w2] = size(WLtilde); dual_wav_props.A_L = WLtilde - x2(1:w1,1:w2);
    [w1, w2] = size(WRtilde); dual_wav_props.A_R = WRtilde - x2((M-w1+1):M,(M-w2+1):M); 
end





function [h0,h1,g0,g1]=compute_spline_filters(N, Ntilde)
  Navg=(N+Ntilde)/2;
  vals=computeQN(Navg);
  
  h0 = 1;
  for k=1:(N/2)
    h0=conv(h0,[1/4 1/2 1/4]);
  end
  h0=h0*vals(1);
  h0 = conv(h0, vals);
  
  g0=1;
  for k=1:(Ntilde/2)
    g0=conv(g0,[1/4 1/2 1/4]);
  end
  
  x = sqrt(2)/abs(sum(h0));
  g0=g0/x;
  h0=h0*x;
  
  h1=g0.*(-1).^((-(length(g0)-1)/2):((length(g0)-1)/2));
  g1=h0.*(-1).^((-(length(h0)-1)/2):((length(h0)-1)/2));
end