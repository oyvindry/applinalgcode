function [W,A_pre_inv]=bw_compute_left_ortho(g0,g1,N,K)
    R = N; L = -N+1;
    
    [C,C_c]=findc(R,K,N); % dimension of C is (2N-1)xN
    Cpinv = inv(C'*C)*C';
    
    Gup = Gsegment(g0,L:R,(-R+1):(K-1),(-2*R+2):2:(2*K-2));
    Gdown = Gsegment(g0,L:R,K:(2*K+R-2),(-2*R+2):2:(2*K-2));
    
    C = double(C);
    Cpinv = double(Cpinv);
    
    % Theorem 3.3
    X_e = Cpinv*Gup*C;
    Z_e = Gdown*C;
    
    % Step 1: Make phi-supports staggered.
    [Qmatr,Rmatr] = qr((flipud(C((R+K-N):end,:)))');
    P = fliplr(Qmatr);
    invP = flipud(Qmatr');

    % Theorem 3.4
    X_e = invP*X_e*P;
    Z_e = Z_e*P;
    A_pre_inv = C((R+K-N):end,:)*P;
    
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
    % A_pre_inv
    A_pre_inv = double(A_pre_inv);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                            Psi-functions                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    val = 3*N - K - 1;
    if mod(val,2) == 1
        val = val + 1;
    end
    N0 = val/2;
    T = N0 - 1;
    S=1:(N0+K-N);
    
    % Step 3: Project the psi-functions onto the orthogonal complement of the,
    % phi-functions. 
    
    newcols = Gsegment(g0,(-N+1):N,0:(2*T + K),(2*N):2:(2*T));
    [Rmatr, newcols] = expand_cols_smallest([X_e;Z_e], newcols);
    G = [Rmatr newcols];
    
    
    lastmatr = G*(G(S,:))';
    idmatr = eye(N0+K-N);
    [idmatr,lastmatr]=expand_cols_smallest(idmatr,lastmatr);
    X_o = idmatr - lastmatr;

    % Step 4: Make psi-supports staggered
    % [Qmatr,Rmatr] = qr((flipud(Z_o(1:2:(2*N-1),:)))');
    % P = fliplr(Qmatr);
    [Rmatr,jb] = rref(double((flipud(X_o))'));
    jb = sort(size(X_o,1) + 1 - jb);
    [Lmatr,Umatr] = lu((flipud(X_o(jb,:)))');
    invP = flipud(Lmatr');
    P = inv(invP);
    X_o = X_o*P;

    % Step 5: Orthogonalize psi-functions    
    Y = X_o'*X_o;
    P = inv(chol(Y));
    X_o = X_o*P;
    
    % Assemble W
    X_e = [X_e; Z_e];
    [X_e, X_o] = expand_cols_smallest(X_e, X_o);
    numcolsW = K + N;
    W = zeros(size(X_e,1),numcolsW);
    W(:,1:(K-N)) = X_o(:,1:(K-N));                   % K-N psi-functions at the beginning.
    W(:,K - N     + 2*(1:N0)) = X_o(:,(K-N+1):end);  % the remaining psi-functions
    W(:,K - N - 1 + 2*(1:N))  = X_e;                 % all phi functions
    
    % Add internal phi-functions
    insertphi = Gsegment(g0, (-N+1):N,  0:(K-N + 2*N0-2 + N), K-N + 2*(N:(N0-1)) );
    [W, insertphi] = expand_cols_smallest(W, insertphi);
    W( :, K-N + 2*(N:(N0-1)) + 1) = insertphi;
        
    % W'*W    
end    

function [C,C_c]=findc(R,K,N)
    C_c = sym(zeros(N));
    C_c(1,1) = sym(1);

    L1 = 2/sym((R+K-1)); L0 = sym((R-K)/(R+K-1));
    C_c(1,2) = L0; C_c(2,2) = L1;

    for n = 1:(N-2) % fill out column n+2. n is as in paper
        betaval = (sym(n^2)/sym(4*n^2-1))*(1-sym(n^2)/sym((R+K-1)^2)) ;
        C_c(1,n+2) = L0*C_c(1,n+1); 
        C_c(2:N,n+2) = L1*C_c(1:(N-1),n+1) + L0*C_c(2:N,n+1);
        C_c(:,n+2) = C_c(:,n+2) - betaval*C_c(:,n);
    end

    C_0 = sym(zeros(R+K-1,N));
    for k=0:(N-1)
        C_0(:,k+1) = sym((-R+1):(K-1)).^k;
    end

    C = C_0*C_c;
    % C'*C check for orthogonality
end

function val=Gsegment(g0,supp,rowrange,colrange)
    val = zeros(length(rowrange),length(colrange));
    k=1;
    for col_ind = colrange
        actualinds =  supp + col_ind;
        [intersec,i1,i2] = intersect(rowrange,actualinds);
        val(i1,k) = g0( actualinds(i2) - actualinds(1) + 1 );
        k = k+1;
    end
end

function [Anew,Bnew]=expand_cols_smallest(A,B)
    if size(A,1) > size(B,1)
        Anew = A;
        Bnew = [ B; zeros(size(A,1)-size(B,1),size(B,2))];
    else
        Anew = [ A; zeros(size(B,1)-size(A,1),size(A,2))];
        Bnew = B;
    end
end