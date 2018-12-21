function [wav_props, dual_wav_props]=find_wav_props(wave_name, m, bd_mode, length_signal)
    % Computes the properties of a wavelet with the given name. What properties 
    % are computed depend on the bd_mode parameter, m, and length_signal.
    %
    % wave_name: Name of the wavelet. Possible names are:
    %            'cdf97' - CDF 9/7 wavelet
    %            'spline53' - Spline 5/3 wavelet
    %            'splinex.x' - Spline wavelet with given number of vanishing 
    %                          moments for each filter
    %            'pwl0'  - Piecewise linear wavelet with 0 vanishing moments
    %            'pwl2'  - Piecewise linear wavelet with 2 vanishing moments
    %            'Haar'  - The Haar wavelet
    %            'dbX'   - Daubechies orthnormal wavelet with X vanishing
    %                      moments
    %            'symX'  - Symmlets: A close to symmetric, orthonormal wavelet 
    %                      with X vanishing moments
    % m:         Number of resolutions. Default: 1
    % bd_mode:   Boundary extension mode. Possible modes are. 
    %            'per'    - Periodic extension
    %            'symm'   - Symmetric extension (default)
    %            'none'   - Take no extra action at the boundaries
    %            'bd'     - Boundary wavelets
    % length_signal: Length of the input signal. Default: 0.

    if (~exist('m','var')) m = 1; end
    if (~exist('bd_mode','var')) bd_mode = 'symm'; end
    if (~exist('length_signal','var')) length_signal = 0; end
    wav_props.wave_name = wave_name; dual_wav_props.wave_name = wave_name;
    wav_props.m = m; dual_wav_props.m = m;
    wav_props.length_signal = length_signal; dual_wav_props.length_signal = length_signal;
    wav_props.offset_L = 0; dual_wav_props.offset_L = 0;
    wav_props.offset_R = 0; dual_wav_props.offset_R = 0;

    if strcmpi(wave_name,'Haar')
        % do nothing
    elseif strcmpi(wave_name(1:2), 'db')
        N = str2double(wave_name(3:end));
        if N > 1 
            [wav_props,dual_wav_props] = wav_props_ortho(N, wav_props, dual_wav_props, bd_mode);
        end
    elseif strcmpi(wave_name(1:3), 'sym')
        N = str2double(wave_name(4:end));
        if N > 1
            [wav_props,dual_wav_props] = wav_props_ortho(N, wav_props, dual_wav_props, bd_mode, 1);
        end
    elseif strcmpi(wave_name, 'pwl0')
        [wav_props, dual_wav_props] = wav_props_pwl0(wav_props, dual_wav_props, bd_mode);
    elseif strcmpi(wave_name, 'pwl2')
        [wav_props, dual_wav_props] = wav_props_pwl2(wav_props, dual_wav_props, bd_mode);
    elseif strcmpi(wave_name, 'spline53')
        [wav_props, dual_wav_props] = wav_props_53(wav_props, dual_wav_props, bd_mode);
    elseif strcmpi(wave_name, 'cdf97')
        [wav_props, dual_wav_props] = wav_props_97(wav_props, dual_wav_props, bd_mode);
    elseif strcmpi(wave_name(1:6), 'spline')
        N = str2double(wave_name(7));
        Ntilde = str2double(wave_name(9));
        [dual_wav_props.g0, dual_wav_props.g1, wav_props.g0, wav_props.g1]=compute_spline_filters(N, Ntilde);
        if strcmpi(bd_mode, 'bd')
            [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde]=wav_props_biortho_bd(N, Ntilde, wav_props, dual_wav_props);
            % swap filters if offset is odd
            if mod(wav_props.offset_L,2) == 1
                g0temp = wav_props.g0; wav_props.g0 = wav_props.g1; wav_props.g1 = g0temp;
                g0temp = dual_wav_props.g0; dual_wav_props.g0 = dual_wav_props.g1; dual_wav_props.g1 = g0temp;
            end
            [wav_props.A_L,wav_props.A_R]=find_AL_AR(WL, WR, wav_props);
            [dual_wav_props.A_L,dual_wav_props.A_R]=find_AL_AR(WLtilde, WRtilde, dual_wav_props);
        end 
    end
end

function [wav_props, dual_wav_props]=wav_props_pwl0(wav_props, dual_wav_props, bd_mode)
    wav_props.last_even = 0;
    wav_props.lambdas = 0.5;
    wav_props.alpha = sqrt(2);
    wav_props.beta = sqrt(2);
    dual_wav_props.g0 = sqrt(2);
    dual_wav_props.g1 = [-1/2 1 -1/2]*sqrt(2);
    wav_props.g0 = [1/2 1 1/2]/sqrt(2);
    wav_props.g1 = 1/sqrt(2);
    
    WL = 0; WLtilde = 0; WR = 0; WRtilde = 0;
    if strcmpi(bd_mode, 'bd')
        [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde] = wav_props_biortho_bd(2, 0, wav_props, dual_wav_props); % TODO: test. how many vanishing moments here?
    end 
    
    [wav_props, dual_wav_props]=set_wav_props(wav_props, dual_wav_props, bd_mode, WL, WLtilde, WR, WRtilde);
end

function [wav_props, dual_wav_props]=wav_props_pwl2(wav_props, dual_wav_props, bd_mode)
    wav_props.last_even = 0;
    wav_props.lambdas = [0.25; 0.5];
    wav_props.alpha = sqrt(2);
    wav_props.beta = sqrt(2);
    dual_wav_props.g0 = [-1/8 1/4 3/4 1/4 -1/8]*sqrt(2);
    dual_wav_props.g1 = [-1/2 1 -1/2]*sqrt(2);
    wav_props.g0 = [1/2 1 1/2]/sqrt(2);
    wav_props.g1 = [-1/8 -1/4 3/4 -1/4 -1/8]/sqrt(2);
    
    WL = 0; WLtilde = 0; WR = 0; WRtilde = 0;
    if strcmpi(bd_mode, 'bd')
        [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde] = wav_props_biortho_bd(2, 2, wav_props, dual_wav_props);
    end 
    
    [wav_props, dual_wav_props]=set_wav_props(wav_props, dual_wav_props, bd_mode, WL, WLtilde, WR, WRtilde);
end


function [wav_props, dual_wav_props]=wav_props_53(wav_props, dual_wav_props, bd_mode)
    wav_props.last_even = 0;
    wav_props.lambdas = [-1; 0.125];
    wav_props.alpha = 2;
    wav_props.beta = 0.5;
    dual_wav_props.g0 = [-1/4 1/2 3/2 1/2 -1/4];
    dual_wav_props.g1 = [-1/4 1/2 -1/4];
    wav_props.g0 = [1/4 1/2 1/4];
    wav_props.g1 = [-1/4 -1/2 3/2 -1/2 -1/4];
    
    WL = 0; WLtilde = 0; WR = 0; WRtilde = 0;
    if strcmpi(bd_mode, 'bd')
        [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde] = wav_props_biortho_bd(2, 2, wav_props, dual_wav_props);
    end 
    
    [wav_props, dual_wav_props]=set_wav_props(wav_props, dual_wav_props, bd_mode, WL, WLtilde, WR, WRtilde);
end

function [wav_props, dual_wav_props]=wav_props_97(wav_props, dual_wav_props, bd_mode)
    wav_props.last_even = 0;
    [wav_props.lambdas, wav_props.alpha, wav_props.beta, dual_wav_props.g0, dual_wav_props.g1, wav_props.g0, wav_props.g1] = lifting_fact_97();

    WL = 0; WLtilde = 0; WR = 0; WRtilde = 0;
    if strcmpi(bd_mode, 'bd')
        [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde] = wav_props_biortho_bd(4, 4, wav_props, dual_wav_props);
    end 
    
    [wav_props, dual_wav_props]=set_wav_props(wav_props, dual_wav_props, bd_mode, WL, WLtilde, WR, WRtilde);
end

function [wav_props, dual_wav_props]=set_wav_props(wav_props, dual_wav_props, bd_mode, WL, WLtilde, WR, WRtilde)
    if strcmpi(bd_mode, 'bd')
        if mod(wav_props.offset_L,2) == 1
            wav_props.last_even = ~wav_props.last_even;
            betaval = wav_props.alpha; wav_props.alpha = wav_props.beta; wav_props.beta = betaval;
        end
    end 
    
    dual_wav_props.lambdas = -wav_props.lambdas; 
    dual_wav_props.alpha = 1/wav_props.alpha; 
    dual_wav_props.beta = 1/wav_props.beta;
    dual_wav_props.last_even = ~wav_props.last_even;
    
    if strcmpi(bd_mode, 'bd')            
        [wav_props.A_L,wav_props.A_R]=find_AL_AR_lifting(WL, WR, wav_props);
        [dual_wav_props.A_L,dual_wav_props.A_R]=find_AL_AR_lifting(WLtilde, WRtilde, dual_wav_props);
    end 
end

function [wav_props, dual_wav_props]=wav_props_ortho(N, wav_props, dual_wav_props, bd_mode, type)
    % N:    Number of vanishing moments
    % type: The type of orthonormal wavelet.
    %       0: Daubechies wavelets with minimum phase (default)  
    %       1: Symmlets - wavelets with close to linear phase (almost symmetric)
    % 
    if (nargin == 4)
        type = 0;
    end
        
    h0 = 0; h1 = 0;
    if (type == 0)
        [h0, h1, wav_props.g0, wav_props.g1] = h0h1_compute_ortho(N);
    elseif (type == 1)
        [h0, h1, wav_props.g0, wav_props.g1] = h0h1_compute_sym(N);
    end
    [wav_props.lambdas, wav_props.alpha, wav_props.beta, wav_props.last_even] = lifting_fact_ortho(h0, h1);
    
    if strcmpi(bd_mode, 'bd')
        [wav_props, WL, WR] = wav_props_ortho_bd(N, wav_props);
        
        if mod(wav_props.offset_L, 2) == 1
            wav_props.last_even = ~wav_props.last_even;
            betaval = wav_props.alpha; wav_props.alpha = wav_props.beta; wav_props.beta = betaval;
        end
    end
        
    dual_wav_props.lambdas = -fliplr(wav_props.lambdas);
    dual_wav_props.alpha = 1/wav_props.alpha;
    dual_wav_props.beta = 1/wav_props.beta;
    dual_wav_props.last_even = ~wav_props.last_even;
    
    if strcmpi(bd_mode, 'bd')            
        [wav_props.A_L, wav_props.A_R]=find_AL_AR_lifting(WL, WR, wav_props);
        [dual_wav_props.A_L,dual_wav_props.A_R]=find_AL_AR_lifting(WL, WR, dual_wav_props);
        dual_wav_props.offset_L = wav_props.offset_L; 
        dual_wav_props.offset_R = wav_props.offset_R; 
        dual_wav_props.A_R_pre_inv = wav_props.A_R_pre_inv;
        dual_wav_props.A_L_pre_inv = wav_props.A_L_pre_inv;
    end 
end

function [wav_props, WL, WR]=wav_props_ortho_bd(N, wav_props)
    % Compute K_L and K_R
    K_L = N; K_R = N;
    s = mod(wav_props.length_signal, 2^wav_props.m);
    if s > 0
        toadd = 2^wav_props.m - s;
        K_L = K_L + floor(toadd/2);
        K_R = K_R + ceil(toadd/2);
    end
    wav_props.offset_L = K_L - N; wav_props.offset_R = K_R - N;
    Mint = (wav_props.length_signal -2*N + K_L + K_R)/2^wav_props.m;
    assert(K_L+K_R+1 <= Mint)
        
    % First the right edge
    wav_props.g0 = flip(wav_props.g0); wav_props.g1 = flip(wav_props.g1);
    [W, A_pre_inv] = bw_compute_left_ortho(wav_props.g0, wav_props.g1, N, K_R);

    % Mirror right hand side variables
    wav_props.A_R_pre_inv = fliplr(flipud(A_pre_inv));
    %A_R_pre_inv = wav_props.A_R_pre_inv;
    WR = fliplr(flipud(W));
    for k = (size(W,2)-(K_R-N)):(-2):1
        WR(:, [k-1 k]) = WR(:, [k k-1]); 
    end



    % Then the left edge
    wav_props.g0 = flip(wav_props.g0); wav_props.g1 = flip(wav_props.g1);
    [WL, wav_props.A_L_pre_inv] = bw_compute_left_ortho(wav_props.g0, wav_props.g1, N, K_L);

    %A_L_pre_inv = wav_props.A_L_pre_inv;
    %h0=wav_props.g0;
    %save(sprintf('DBFilter/h0_db%d.mat', N), 'h0');

end

function [A_L,A_R]=find_AL_AR_lifting(WL, WR, wav_props)
    M = max(size(WL)+size(WR));
    if mod(M-wav_props.length_signal,2) == 1
        M = M + 1;
    end
    
    x1 = eye(M);
    [f, prefilter]=find_kernel(wav_props, 0, 0, 0, 0, 'none');
    x1 = idwt1_impl_internal(x1, f, 1, 'none', prefilter, [wav_props.offset_L wav_props.offset_R], 'time');
    [w1, w2] = size(WL); A_L = WL - x1(1:w1,1:w2);
    %matr = zeros(size(x1)); 
    %matr(1:w1,1:w2) = A_L;
    [w1, w2] = size(WR); A_R = WR - x1((M-w1+1):M,(M-w2+1):M);
    % A_L
    % A_R
    A_L = double(A_L); A_R = double(A_R);
    %matr((M-w1+1):M,(M-w2+1):M) = A_R;
    %x1+matr
    %(x1+matr)'*(x1+matr)
    %x1
    %matr
end

function [wav_props, dual_wav_props, WL, WLtilde, WR, WRtilde]=wav_props_biortho_bd(N, Ntilde, wav_props, dual_wav_props)
    Nprime = max(N, Ntilde);
    R = (length(wav_props.g0)-1)/2; Rtilde = (length(dual_wav_props.g0)-1)/2;
    L  = -R; Ltilde = -Rtilde;
        
    K_L = -L; K_L_tilde = K_L + Ntilde - N;
    if K_L_tilde < -Ltilde
        K_L_tilde = -Ltilde;
        K_L = K_L_tilde + N - Ntilde;
    end
    K_R = K_L; K_R_tilde = K_L_tilde;
        
    s = wav_props.length_signal + R + L - 2*Nprime + K_L + K_R - 1;
    if mod(s,2^wav_props.m) > 0
        toadd = 2^wav_props.m - mod(s,2^wav_props.m);
        K_L = K_L + floor(toadd/2); K_L_tilde = K_L_tilde + floor(toadd/2);
        K_R = K_R + ceil(toadd/2);  K_R_tilde = K_R_tilde + ceil(toadd/2);
    end
    wav_props.offset_L = K_L - N; wav_props.offset_R = K_R - N; dual_wav_props.offset_L = wav_props.offset_L; dual_wav_props.offset_R = wav_props.offset_R;
    Mint = (wav_props.length_signal -2*Nprime + K_L + K_R-1)/2^wav_props.m;
    assert(K_L+K_R <= Mint)
        
        
    [WL, WLtilde, wav_props.A_L_pre_inv, dual_wav_props.A_L_pre_inv] = bw_compute_left_biortho(wav_props.g0, wav_props.g1, N, dual_wav_props.g0, dual_wav_props.g1, Ntilde, K_L, K_L_tilde);
    [WR, WRtilde, wav_props.A_R_pre_inv, dual_wav_props.A_R_pre_inv] = bw_compute_left_biortho(wav_props.g0, wav_props.g1, N, dual_wav_props.g0, dual_wav_props.g1, Ntilde, K_R, K_R_tilde);
    
    % Mirror right hand side variables
    wav_props.A_R_pre_inv = flipud(fliplr(wav_props.A_R_pre_inv));
    dual_wav_props.A_R_pre_inv = flipud(fliplr(dual_wav_props.A_R_pre_inv));
    
    WR = flipud(fliplr(WR)); WRtilde = flipud(fliplr(WRtilde));
end

function [A_L,A_R]=find_AL_AR(WL, WR, wav_props)
    R = (length(wav_props.g0)-1)/2; Rtilde = (length(wav_props.g1)-1)/2;
    L  = -R; Ltilde = -Rtilde;
    M = max(size(WL)+size(WR));
    if mod(M-wav_props.length_signal,2) == 1
        M = M + 1;
    end
    x1 = zeros(M);
    x1(:,1:2:end)=Gsegment(wav_props.g0, L:R, 0:(M-1), 0:2:(M-1));
    x1(:,2:2:end)=Gsegment(wav_props.g1, Ltilde:Rtilde, 0:(M-1), 1:2:(M-1));
    
    [w1, w2] = size(WL); A_L = WL - x1(1:w1,1:w2);
    [w1, w2] = size(WR); A_R = WR - x1((M-w1+1):M,(M-w2+1):M);
    % A_L
    % A_R
    A_L = double(A_L); A_R = double(A_R);
end



function [h0,h1,g0,g1]=compute_spline_filters(N, Ntilde)
  Navg=(N+Ntilde)/2;
  vals=compute_QN(Navg);
  
  h0 = 1;
  for k=1:(N/2)
    h0=conv(h0,[1/4 1/2 1/4]);
  end
  %h0=h0*vals(1);
  h0 = conv(h0, vals);
  
  g0=1;
  for k=1:(Ntilde/2)
    g0=conv(g0,[1/4 1/2 1/4]);
  end
  
  %x = sqrt(2)/abs(sum(h0));
  %g0=g0/x;
  %h0=h0*x;
  
  h1=g0.*(-1).^((-(length(g0)-1)/2):((length(g0)-1)/2));
  g1=h0.*(-1).^((-(length(h0)-1)/2):((length(h0)-1)/2));
end

function val=Gsegment(g0,supp,rowrange,colrange)
    val = sym(zeros(length(rowrange),length(colrange)));
    k=1;
    for col_ind = colrange
        actualinds =  supp + col_ind;
        [intersec,i1,i2] = intersect(rowrange,actualinds);
        val(i1,k) = g0( actualinds(i2) - actualinds(1) + 1 );
        k = k+1;
    end
end

function [lambdas, alpha, beta, h0, h1, g0, g1]=lifting_fact_97()
    [h0,h1,g0,g1] = h0h1_compute_97(); % Should have 9 and 7 filter coefficients.
    h00 = h0(1:2:9);
    h01 = h0(2:2:9);
    h10 = h1(1:2:7);
    h11 = h1(2:2:7); 
    
    lambdas=zeros(4,1);
    
    lambdas(1) = -h00(1)/h10(1);
    h00(1:5) = h00(1:5) + conv(h10(1:4),[lambdas(1) lambdas(1)]);
    h01(1:4) = h01(1:4) + conv(h11(1:3),[lambdas(1) lambdas(1)]); 
    
    lambdas(2) = -h10(1)/h00(2);
    h10(1:4) = h10(1:4)+conv(h00(2:4),[lambdas(2) lambdas(2)]);
    h11(1:3) = h11(1:3)+conv(h01(2:3),[lambdas(2) lambdas(2)]);
    
    lambdas(3) = -h00(2)/h10(2);
    h00(2:4) = h00(2:4)+conv(h10(2:3),[lambdas(3) lambdas(3)]);
    h01(2:3) = h01(2:3)+conv(h11(2:2),[lambdas(3) lambdas(3)]); 
    
    lambdas(4) = -h10(2)/h00(3);
    h10(2:3) = h10(2:3)+conv(h00(3:3),[lambdas(4) lambdas(4)]);
    
    alpha = h00(3);
    beta  = h11(2);
end

function [h0,h1,g0,g1]=h0h1_compute_97()
  N=4;
  vals=compute_QN(N);
   
  rts=roots(vals)';
  rts1=rts(find(abs(imag(rts))>0.001)); % imaginary roots
  rts2=rts(find(abs(imag(rts))<0.001)); % real roots
  
  h0=1;
  for rt=rts1
    h0=conv(h0,[-rt 1]);
  end
  for k=1:(N/2)
    h0=conv(h0,[1/4 1/2 1/4]);
  end
  h0=h0*vals(1);
  
  g0=1;
  for rt=rts2
    g0=conv(g0,[-rt 1]);
  end
  for k=1:(N/2)
    g0=conv(g0,[1/4 1/2 1/4]);
  end
  
  
  h0=real(h0);
  g0=real(g0);  
  x = sqrt(2)/abs(sum(h0));
  g0=g0/x;
  h0=h0*x;
  
  h1=g0.*(-1).^((-(length(g0)-1)/2):((length(g0)-1)/2));
  g1=h0.*(-1).^((-(length(h0)-1)/2):((length(h0)-1)/2));
end

function vals=compute_QN(N)
    % Compute the coefficients in Q^(N)((1-cos(w))/2)
    k=0:(N-1);
    QN = 2*factorial(N+k-1)./(factorial(k).*factorial(N-1));
    vals=zeros(1,2*N-1);
    vals=QN(1);
    start=1;
    for k=2:N
        start=conv(start,[-1/4 1/2 -1/4]);
        vals=[0 vals 0]+QN(k)*start;
    end
end

function [h0, h1, g0, g1]=h0h1_compute_ortho(N)
    % Comptues the wavelet coefficients of the orthonormal Daubechies wavelet
    % N vanishing moments and with minimum phase   
    vals=compute_QN(N);
    rts=roots(vals)';
    rts1=rts(find(abs(rts)>1));

    g0=1;
    for rt=rts1
        g0=conv(g0,[-rt 1]);
    end
    g0 = real(g0);
    K=sqrt(vals(1)*(-1)^(length(rts1))/abs(prod(rts1)));
    g0=K*g0;
    for k=1:N
        g0=conv(g0,[1/2 1/2]);
    end
    
    % Ensuring integral is positive - This part of the code requiere some more
    % testing
    if (sum(g0) < 0)
        g0 = -g0;
    end
    h0=fliplr(g0);
    g1=h0.*(-1).^(0:(length(g0)-1)); % It seems to me that this should be 
                                     % multiplied by -1 
    h1=fliplr(g1);
end

function [h0, h1, g0, g1]=h0h1_compute_sym(N)
    % Comptues the wavelet coefficients of the orthonormal wavelet with N
    % vanishing moments and close to linear phase. This makes the wavelet
    % almost symmetric. These wavelets are called 'symlets'
    %
    % This function relies on matlabs wavelet coefficients. In the next version
    % this will be changed 
    
    currDWTmode = dwtmode('status', 'nodisp');
    dwtmode('per','nodisp');
    nu = 7;
    n = 2^nu;
    x = zeros([1,n]);
    x(ceil(N/2)) = 1;
    
    S = [2^(nu-1); 2^(nu-1); n]; % compute the S given by wavedec
    wave_name = sprintf('sym%d', N);
    
    y = waverec(x, S, wave_name);
    if (mod(N,2) == 1) % is odd
        g0 = y(1:2*N);
    else % is even 
        g0 = [y(end), y(1:2*N-1)];
    end
    
    h0=fliplr(g0);
    g1=h0.*(-1).^(0:(length(g0)-1)); 
    h1=fliplr(g1);
    
    dwtmode(currDWTmode, 'nodisp');
end

function [lambdas, alpha, beta, last_even] = lifting_fact_ortho(h0, h1)
    % Compute a lifting factorization so that
    % [alpha 0; beta 0] = \Lambda_n \cdots \Lambda_1 H
    % lambdas: [\Lambda_1; \Lambda_2; \cdots \Lambda_n]
    % last_even: If \Lambda_n is even
    stepnr=1;
    len1=length(h0)/2; len2=len1;
    lambdas=zeros(len1+1,2);
    if mod(len1,2)==0
        h00=h0(1:2:length(h0));
        h01=h0(2:2:length(h0));
        h10=h1(1:2:length(h1));
        h11=h1(2:2:length(h1));
  
        lambda1=-h00(1)/h10(1);
        h00=h00+lambda1*h10; 
        h01=h01+lambda1*h11;
        start1=2; end1=len1; len1=len1-1; start2=1; end2=len2;
        lambdas(stepnr,:)=[lambda1 0];
    else
        h00=h0(2:2:length(h0));
        h01=h0(1:2:length(h0));
        h10=h1(2:2:length(h1));
        h11=h1(1:2:length(h1));
    
        lambda1=-h10(len1)/h00(len1); 
        h10=h10+lambda1*h00; 
        h11=h11+lambda1*h01;
        start2=1; end2=len2-1; len2=len2-1; start1=1; end1=len1;
        lambdas(stepnr,:)=[0 lambda1];
    end
  
    %[h00 h01; h10 h11]
    %conv(h00,h11)-conv(h10,h01)
    stepnr=stepnr+1;

    %[h00 h01; h10 h11]
    %conv(h00,h11)-conv(h10,h01)
    while len2>0 % Stop when the second element in the first column is zero
        if len1>len2 % Reduce the degree in the first row. 
            lambda1=-h00(start1)/h10(start2);
            lambda2=-h00(end1)/h10(end2);
            h00(start1:end1)=h00(start1:end1)+conv(h10(start2:end2),[lambda1 lambda2]);
            h01(start1:end1)=h01(start1:end1)+conv(h11(start2:end2),[lambda1 lambda2]);
            start1=start1+1; end1=end1-1; len1=len1-2;
        else % reduce the degree in the second row. 
            lambda1=-h10(start2)/h00(start1);
            lambda2=-h10(end2)/h00(end1);
            h10(start2:end2)=h10(start2:end2)+conv(h00(start1:end1),[lambda1 lambda2]);
            h11(start2:end2)=h11(start2:end2)+conv(h01(start1:end1),[lambda1 lambda2]);
            start2=start2+1; end2=end2-1; len2=len2-2;
        end
        lambdas(stepnr,:)=[lambda1 lambda2];
        stepnr=stepnr+1;
    
        %[h00 h01; h10 h11]
        %conv(h00,h11)-conv(h10,h01)
    end
  
    % Add the final lifting, and alpha,beta
    alpha=sum(h00);
    beta=sum(h11);
    lastlift=-sum(h01)/beta;
    if mod(length(h0)/2,2)==0
        lambdas(stepnr,:)=[0 lastlift];
    else
        lambdas(stepnr,:)=[lastlift 0];
    end
    last_even = 1;
    %[h00 h01; h10 h11]
end
