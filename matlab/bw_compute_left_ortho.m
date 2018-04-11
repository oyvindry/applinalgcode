function [W,A_pre,A_pre_inv]=bw_compute_left_ortho(g0,N,debugMode)
    if (nargin < 3) 
        debugMode = 0;
    end
    R = N;
    
    [C,C_c]=findc(R,N); % dimension of C is (2N-1)xN
    Cpinv = inv(C'*C)*C';
    
    Gup =   Geven(g0,(-N+1):N,(-N+1):(KN-1),(-2*N+2):2:(2*N-2));
    Gdown = Geven(g0,(-N+1):N,N:(3*N-2)    ,(-2*N+2):2:(2*N-2));
    
    C = double(C);
    Cpinv = double(Cpinv);
    
    % Theorem 3.3
    X_e = Cpinv*Gup*C; % dimensions NxN
    Z_e = Gdown*C; % dimensions (N+R-1)xN
    
    % Step 1: Make phi-supports staggered.
    [Q,R] = qr((flipud(C(N:end,:)))');
    P = fliplr(Q);
    invP = flipud(Q');

    % Theorem 3.4
    X_e = invP*X_e*P;
    Z_e = Z_e*P;
    A_pre_inv = C(R:end,:)*P;
    
    % Step 2: Orthogonalize phi-functions
    ls = eye(N^2)-kron(X_e',X_e');
    rs = reshape(Z_e'*Z_e, [N^2, 1]); % From matrix to vector
    Y  = reshape(ls\rs, [N, N]);      % From vector to matrix
    %fprintf('cond(Y): %g\n', cond(Y))

    P = inv(chol(Y));
    
    % Theorem 3.4
    X_e = inv(P)*X_e*P;
    Z_e = Z_e*P;
    
    A_pre_inv = A_pre_inv*P;
    A_pre = inv(A_pre_inv);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                            Psi-functions                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    P = inv(chol(Y));
    
    % Theorem 4.1
    X_o = X_o*P;
    Z_o = Z_o*P;
    
    W = zeros(3*N-1,2*N);
    W(:, 1:2:(2*N-1)) = [X_e; Z_e];
    W(:, 2:2:(2*N)) = [X_o; Z_o];
end    

function [C,C_c]=findc(R,N)
    C_c = sym(zeros(N));
    C_c(1,1) = sym(1);

    L1 = 2/sym((R+N-1)); L0 = sym((R-N)/(R+N-1));
    C_c(1,2) = L0; C_c(2,2) = L1;

    for n = 1:(N-2) % fill out column n+2. n is as in paper
        betaval = (sym(n^2)/sym(4*n^2-1))*(1-sym(n^2)/sym((R+N-1)^2)) ;
        C_c(1,n+2) = L0*C_c(1,n+1); 
        C_c(2:N,n+2) = L1*C_c(1:(N-1),n+1) + L0*C_c(2:N,n+1);
        C_c(:,n+2) = C_c(:,n+2) - betaval*C_c(:,n);
    end

    C_0 = sym(zeros(R+N-1,N));
    for k=0:(N-1)
        C_0(:,k+1) = sym((-R+1):(N-1)).^k;
    end

    C = C_0*C_c;
    % C'*C check for orthogonality
end

function val=Geven(g0,supp,rowrange,colrange)
    val = sym(zeros(length(rowrange),length(colrange)));
    k=1;
    for col_ind = colrange
        actualinds =  supp + col_ind;
        [intersec,i1,i2] = intersect(rowrange,actualinds);
        val(i1,k) = g0( actualinds(i2) - actualinds(1) + 1 );
        k = k+1;
    end
end