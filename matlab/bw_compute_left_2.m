function [W,A_pre,A_pre_inv]=bw_compute_left_2()

    global wavlib_h0 wavlib_h1 wavlib_g0 wavlib_g1
    N = length(wavlib_h0)/2;
    % We use the symbolic toolbox for computations with rational numbers 
    C_c = sym(zeros(N));
    C_c(1,1) = sym(1);

    L = 2/sym((2*N-1));
    C_c(2,2) = L;

    for n = 1:(N-2) % fill out column n+2. n is as in paper
        betaval = (sym(n^2)/sym(4*n^2-1))*(1-sym(n^2)/sym((2*N-1)^2)) ;
        C_c(2:N,n+2) = L*C_c(1:(N-1),n+1);
        C_c(:,n+2) = C_c(:,n+2) - betaval*C_c(:,n);
    end

    C_0 = sym(zeros(2*N-1,N));
    for k=0:(N-1)
        C_0(:,k+1) = sym((-N+1):(N-1)).^k;
    end

    C = C_0*C_c;
    m = C'*C;
    for k=1:N
        C(:,k) = C(:,k)/sym(sqrt(m(k,k)));
    end

    % dimension of C is (2*N-1)xN
    C = double(C); % Use numerically from now on
    %C'*C % Check for orthogonality

    % For the left edge we need the ((-N+1):(3N-2))x((-2N+2):(2N-2)) segment of
    % G. The first 2N-1 rows are for X, the remaining 2N-1 rows are for Z.
    % This equals the ((N-1):(5N-4))x(0:(4N-3)) segment of G, due to the
    % Toeplitz structure.
    
    seg1 = zeros(2*N-1);
    seg2 = zeros(2*N-1);
    for i =1:N-1
        seg1(1:2*i, i) = wavlib_g0((end-2*i+1):end); 
    end
    for i = 0:N-1
        seg1((1+2*i):end, N+i) = wavlib_g0(1:(end - 2*i - 1));
        seg2(1:(2*i+1), N+i) = wavlib_g0((end - 2*i):end);
    end
        
    % Step 1: Make phi-supports staggered.
    [Q,R] = qr((flipud(C(N:end,:)))');
    P = fliplr(Q);
    C = C*P;

    ls2 = eye((2*N-1)^2) - kron(seg1',seg1');
    rs2 = reshape(seg2'*seg2, [(2*N-1)^2, 1]);
    Y  = reshape(ls2\rs2, [2*N-1, 2*N-1]);
    Y = C'*Y*C; % Gramm matrix
    
    % Step 2: Orthogonalize phi-functions
    P = ortho_from_gramm(Y, N);
    C = C*P;
    A_pre_inv = C(N:end,:);
    A_pre = inv(A_pre_inv);
    
    % Theorem 3.3
    X_e = pinv(C)*seg1*C; % dimensions NxN
    Z_e =    seg2*C; % dimensions (2N-1)xN
    
    % Step 3: Project the psi-functions onto the orthogonal complement of the,
    % phi-functions. 
    X_o = eye(N) - X_e*X_e';
    Z_o = -Z_e*X_e';

    % Step 4: Make psi-supports staggered
    [Q,R] = qr((flipud(Z_o(1:2:(2*N-1),:)))');
    P = fliplr(Q);
    
    % Theorem 4.1
    X_o = X_o*P;
    Z_o = Z_o*P;

    % Step 5: Orthogonalize psi-functions
    Y = X_o'*X_o + Z_o'*Z_o; % Gram matrix
    P=ortho_from_gramm(Y, N);
    
    % Theorem 4.1
    X_o = X_o*P;
    Z_o = Z_o*P;
    
    W = zeros(3*N-1,2*N);
    W(:, 1:2:(2*N-1)) = [X_e; Z_e];
    W(:, 2:2:(2*N)) = [X_o; Z_o];
end    
    

function P=ortho_from_gramm(Y, N)
    xmatr=zeros(N);
    for k=2:N
        xmatr(1:(k-1),k) = Y(1:(k-1),1:(k-1))\Y(1:(k-1),k);
    end
    P = eye(N) - xmatr;             % The g_n are now orthogonal
    P = P*diag(sqrt(diag(P'*Y*P)).^(-1)); % The g_n are now orthonormal
end