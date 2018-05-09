function [W,A_pre,A_pre_inv, outMat]=bw_compute_left(g0,g1,N,h0,h1,Nt,debugMode)
    if (nargin < 7) 
        debugMode = 0;
    end
    outMat = 0;
    transl = 0;
    if nargin < 4 % orthonormal mode
        R = N;
    else % biorthogonal mode. Translation is needed when the filters don't have dissimilar lengths.
        L = (length(g0)-1)/2;
        Lt = (length(h0)-1)/2;
        Nmin = max(-N,-Nt);
        transl = max(L,Lt) - Nmin;
        R = L+transl;
        Rt = Lt + transl;
    end
    
    [C,C_c]=findc(R,N); % dimension of C is (R+N-1)xN
    Cpinv = inv(C'*C)*C';
    
    if nargin < 4 % orthonormal mode
        [Xgseg,Zgseg]=findsegments(g0,N);
        C = double(C);
        Cpinv = double(Cpinv);
    else
        [Xgseg,Zgseg]=findsegments(g0,N,h0,Nt,0);
    end
    
    % Theorem 3.3
    X_e = Cpinv*Xgseg*C; % dimensions NxN
    Z_e =    Zgseg*C; % dimensions (N+R-1)xN
    % Lemma 4.2
    if nargin >= 4 && Nt > N % biorthogonal mode and need to add more boundary functions
        [X_e,Z_e]=expandxeze(X_e, Z_e, g0, N, h0, Nt)
    end
    outMat = C;
    
    % Step 1: Make phi-supports staggered.
    if nargin < 4 % orthonormal mode
        [Qmatr,Rmatr] = qr((flipud(C(N:end,:)))');
        P = fliplr(Qmatr);
        invP = flipud(Qmatr');
    else % biorthogonal mode. Use LU factorization to obtain exact results
        [Lmatr,Umatr] = lu((flipud(C(R:end,:)))');
        invP = flipud(Lmatr');
        P = inv(invP);
    end

    % Theorem 3.4
    X_e = invP*X_e*P;
    Z_e = Z_e*P;
    A_pre_inv = C(R:end,:)*P;
    
    % Step 2: Orthogonalize phi-functions
    ls = eye(N^2)-kron(X_e',X_e');
    rs = reshape(Z_e'*Z_e, [N^2, 1]); % From matrix to vector
    Y  = reshape(ls\rs, [N, N]);      % From vector to matrix
    Y;
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
    [Qmatr,Rmatr] = qr((flipud(Z_o(1:2:(2*N-1),:)))');
    P = fliplr(Qmatr);
    
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

function [Xgseg,Zgseg]=findsegments(g0,N,h0,Nt,increased)
    R = N; transl = 0; supp = (-N+1):N;
    if nargin >= 3 % biorthogonal mode
        L = (length(g0)-1)/2;
        Lt = (length(h0)-1)/2;
        transl = max(L,Lt)-min(N,Nt)
        R = L+transl;
        supp = ((-L):L)  + transl;
        if increased
            % The columns are (2*N):2:(2*Nt-2) in both cases below
            Xgseg = Geven(g0,supp,N:(Nt-1),(2*N):2:(2*Nt-2); % For X we need the rows N:(Nt-1) of G
            Zgseq = Geven(g0,supp,Nt:(supp(end) + (2*Nt-2)),(2*N):2:(2*Nt-2); % For Z we need the rows >= Nt of G. These end at supp(end) + (2*Nt-2)
            return
        end
    end
    % The columns are (-2*R+2):2:(2*N-2) in both cases below
    Xgseg = Geven(g0,supp,(-R+1):(N-1),(-2*R+2):2:(2*N-2)); % For X we need the rows (-R+1):(N-1) of G
    Zgseq = Geven(g0,supp,N:(2*N+R-2),(-2*R+2):2:(2*N-2)); % For Z we need the rows N:(2*N+R-2) of G
end

function [X_e_new,Z_e_new]=expandxeze(X_e, Z_e, g0, N, h0, Nt)
    [Xgseg,Zgseg]=findsegments(g0, N, h0, Nt, 1);
    X_e_new = [X_e zeros(N,Nt-N); Z_e(1:(Nt-N),:) Xgseg];
        
    toadd = size(Zgseg,1)-size(Z_e((Nt-N+1):end,:),1);
    Z_e_new = [Z_e((Nt-N+1):end,:); zeros(toadd,size(Ze,2))];
    Z_e_new = [Z_e_new Zgseg];
end

function val=Geven(g0,supp,rowrange,colrange)
    val = zeros(length(rowrange),length(colrange));
    k=1;
    for col_ind = colrange
        actualinds =  supp + col_ind;
        [intersec,i1,i2] = intersect(rowrange,actualinds);
        val(i1,k) = g0( actualinds(i2) - actualinds(1) + 1 );
        k = k+1;
    end
end

function val=Godd(g1,supp,rowrange,colrange)
    val = zeros(length(rowrange),length(colrange));
    k=1;
    for col_ind = colrange
        actualinds =  supp + col_ind;
        [intersec,i1,i2] = intersect(rowrange,actualinds);
        val(i1,k) = g1( actualinds(i2) - actualinds(1) + 1 );
        k = k+1;
    end
end