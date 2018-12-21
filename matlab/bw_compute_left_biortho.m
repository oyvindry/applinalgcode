function [W,Wtilde,A_pre_inv,Atilde_pre_inv]=bw_compute_left_biortho(g0, g1, N, h0, h1, Ntilde, K, Ktilde)
    % bw_compute_left_biortho(sym([1/4; 1/2; 1/4]), 2, sym([-1/4; 1/2; 3/2; 1/2; -1/4]), 2); % not correct
    % general: 
    %   [h0,h1,g0,g1]=compute_spline_filters(N, Ntilde); % internal function
    %   bw_compute_left_biortho(g0,g1,N,h0,h1,Ntilde);
    
    A_pre_inv = eye(N);
    Atilde_pre_inv = eye(Ntilde);
    Nprime = max(N,Ntilde);
    R = (length(g0)-1)/2; Rtilde = (length(h0)-1)/2;
    L  = -R; Ltilde = -Rtilde;
    
    val = Nprime + N - K - 1 + max(R,Rtilde);
    if mod(val,2) == 1
        val = val + 1;
    end
    N0 = val/2;
    
    S=[];
    k = 0;
    for indinS=0:(N0 + K - N - 1)
        if k <  K - N + R
           S = [S; k];
           k = k + 1;
        else
           S = [S; k+1];
           k = k + 2;
        end
    end
    S = S + 1;
    
    Stilde=[];
    k = 0;
    for indinS=0:(N0 + K - N - 1)
        if k <  K - N + Rtilde
           Stilde = [Stilde; k];
           k = k + 1;
        else
           S = [Stilde; k+1];
           k = k + 2;
        end
    end
    Stilde = Stilde + 1;
    
    val = 2*N0 - R - Ltilde - 1;
    if mod(val,2) == 1
         val = val + 1;
    end
    T = val/2;
    
    val = 2*N0 - Rtilde - L - 1;
    if mod(val,2) == 1
         val = val + 1;
    end
    Ttilde = val/2;

    % Find initial boundary functions
    Gup = Gsegment(g0,L:R,(-R+1):(K-1),(-2*R+2):2:(2*K-2));
    Gdown = Gsegment(g0,L:R,K:(2*K+R-2),(-2*R+2):2:(2*K-2));
    
    Gtildeup = Gsegment(h0,Ltilde:Rtilde,(-Rtilde+1):(Ktilde-1),(-2*Rtilde+2):2:(2*Ktilde-2));
    Gtildedown = Gsegment(h0,Ltilde:Rtilde,Ktilde:(2*Ktilde+Rtilde-2),(-2*Rtilde+2):2:(2*Ktilde-2));
    
    [C,C_c]=findc(R,K,N);
    [Ctilde,C_ctilde]=findc(Rtilde,Ktilde,Ntilde);
    
    Xtildee = inv(Ctilde'*Ctilde)*Ctilde'*Gtildeup *Ctilde;
    Ztildee = Gtildedown*Ctilde;
    Xe = inv(C'*C)*C'* Gup *C;
    Ze = Gdown*C;
    
    % Make phi staggered
    [Lmatr,Umatr] = lu((flipud(C((R+K-N):end,:)))');
    invP = flipud(Lmatr');
    P = inv(invP);
    Xe = inv(P)*Xe*P;
    Ze = Ze*P;
    if N == Ntilde
        A_pre_inv = C((K+R-N):end,:)*P;
    end
    
    % Make phitilde staggered
    [Lmatr,Umatr] = lu((flipud(Ctilde((Rtilde+Ktilde-Ntilde):end,:)))');
    invP = flipud(Lmatr');
    P = inv(invP);
    Xtildee = inv(P)*Xtildee*P;
    Ztildee = Ztildee*P;
    if N == Ntilde
        Atilde_pre_inv = Ctilde((Ktilde+Rtilde-Ntilde):end,:)*P;
    end
    
    % Address when Ntilde and N are different
    if Ntilde > N    
        XG = Gsegment(g0,L:R,(2*K+L):(Ktilde-1),(2*K):2:(2*Ktilde-2));
        ZG = Gsegment(g0,L:R,Ktilde:(2*Ktilde+R-2),(2*K):2:(2*Ktilde-2));
        
        Xe = [Xe; Ze(1:(Ntilde-N),:)];
        Xe = [Xe [zeros(K+L+N,Ntilde-N); XG]];
        [Ze,ZG] = expand_cols_smallest(Ze((Ntilde-N+1):end,:), ZG);
        Ze = [Ze ZG];
    elseif N > Ntilde
        XG = Gsegment(h0,Ltilde:Rtilde,(2*Ktilde+Ltilde):(K-1),(2*Ktilde):2:(2*K-2));
        ZG = Gsegment(h0,Ltilde:Rtilde,K:(2*K+Rtilde-2),(2*Ktilde):2:(2*K-2));
        
        Xtildee = [Xtildee; Ztildee(1:(N-Ntilde),:)];
        Xtildee = [Xtildee [zeros(Ktilde+Ltilde+Ntilde,N-Ntilde); XG]];
        [Ztildee,ZG] = expand_cols_smallest(Ztildee((N-Ntilde+1):end,:), ZG);
        Ztildee = [Ztildee ZG];
    end
    
    % Bi-orthogonalize phi and phitilde
    ls = eye(Nprime^2) - kron(Xe',Xtildee');
    
    [Ze,Ztildee]=expand_cols_smallest(Ze,Ztildee);
    rs = reshape(Ztildee'*Ze, [Nprime^2, 1]);
    Y  = reshape(ls\rs, [Nprime, Nprime]);
    Y=Y';
    
    [Lmatr,Umatr]= lu(Y);
    P1 = (inv(Lmatr))';
    P2 = inv(Umatr);
    Xe = inv(P1)*Xe*P1;
    Ze = Ze*P1;
    Xtildee = inv(P2)*Xtildee*P2;
    Ztildee = Ztildee*P2;
    
    [Xe; Ze];
    [Xtildee;Ztildee];
    % [Xe' Ze'] * [Xtildee; Ztildee] % ok
    
    if N == Ntilde
        A_pre_inv = A_pre_inv*P1;
        Atilde_pre_inv = Atilde_pre_inv*P2;
        % A_pre_inv
        % Atilde_pre_inv
        A_pre_inv = double(A_pre_inv);
        Atilde_pre_inv = double(Atilde_pre_inv);
    end
    
    % psi-funksjoner

    % Construct Xo
    newcols = Gsegment(g0,L:R,0:(2*T + R + K - N),(2*Nprime):2:(2*T));
    [Rmatr, newcols] = expand_cols_smallest([Xe;Ze], newcols);
    G = [Rmatr newcols];
    
    newcols = Gsegment(h0,Ltilde:Rtilde,0:(2*T + Rtilde + K - N),(2*Nprime):2:(2*T));
    [Rmatr, newcols] = expand_cols_smallest([Xtildee;Ztildee], newcols);
    Gtilde = [Rmatr newcols];
    
    lastmatr = G*(Gtilde(S,:))';
    idmatr = eye(S(N0+K-N));
    idmatr = idmatr(:,S);
    [idmatr,lastmatr]=expand_cols_smallest(idmatr,lastmatr);
    Xo = idmatr - lastmatr;
    
    % Make psi staggered
    [Rmatr,jb] = rref(double((flipud(Xo))'));
    jb = sort(size(Xo,1) + 1 - jb);
    [Lmatr,Umatr] = lu((flipud(Xo(jb,:)))');
    invP = flipud(Lmatr');
    P = inv(invP);
    Xo = Xo*P;
    
    % Construct Xtildeo
    newcols = Gsegment(h0,Ltilde:Rtilde,0:(2*Ttilde + Rtilde + K - N),(2*Nprime):2:(2*Ttilde));
    [Rmatr, newcols] = expand_cols_smallest([Xtildee;Ztildee], newcols);
    Gtilde = [Rmatr newcols];
    
    newcols = Gsegment(g0,L:R,0:(2*Ttilde + R + K - N),(2*Nprime):2:(2*Ttilde));
    [Rmatr, newcols] = expand_cols_smallest([Xe;Ze], newcols);
    G = [Rmatr newcols];
    
    lastmatr = Gtilde*(G(Stilde,:))';
    idmatr = eye(Stilde(N0+K-N));
    idmatr = idmatr(:,Stilde);
    [idmatr,lastmatr]=expand_cols_smallest(idmatr,lastmatr);
    Xtildeo = idmatr - lastmatr;
    
    % Make psitilde staggered
    [Rmatr,jb] = rref(double((flipud(Xtildeo))'));
    jb = sort(size(Xtildeo,1) + 1 - jb);
    [Lmatr,Umatr] = lu((flipud(Xtildeo(jb,:)))');
    invP = flipud(Lmatr');
    P = inv(invP);
    Xtildeo = Xtildeo*P;
    
    % Bi-orthogonalize psi and psitilde
    [Xo,Xtildeo]=expand_cols_smallest(Xo,Xtildeo);
    [Lmatr,Umatr]=lu(Xo'*Xtildeo);
    P1 = (inv(Lmatr))';
    P2 = inv(Umatr);
    Xo = Xo*P1;
    Xtildeo = Xtildeo*P2;

    % Assemble W and Wtilde
    Xe = [Xe; Ze]; Xtildee = [Xtildee; Ztildee];
    W=     assembleW(Xe,      Xo,      g0, L:R,           g1, (-Rtilde):(-Ltilde), K - N, Nprime, N0);
    Wtilde=assembleW(Xtildee, Xtildeo, h0, Ltilde:Rtilde, h1, (-R):(-L),           K - N, Nprime, N0);
    [W, Wtilde]=expand_cols_smallest(W, Wtilde);    
    % only dyadic fractions?
    % Wtilde % only dyadic fractions?
    % W'*Wtilde % Should be identity matrix
end

function W=assembleW(Xe, Xo, g0, suppg0, g1, suppg1, KminN, Nprime, N0)
    [Xe, Xo] = expand_cols_smallest(Xe, Xo);
    numcols = KminN + 2*max(Nprime,N0);
    W = sym(zeros(size(Xe,1),numcols));
    W(:,1:(KminN)) = Xo(:,1:(KminN));                       % K-N psi-functions at the beginning.
    W(:,KminN     + 2*(1:N0))     = Xo(:,(KminN+1):end);  % the remaining psi-functions
    W(:,KminN - 1 + 2*(1:Nprime)) = Xe;                 % all phi functions
    
    if Nprime > N0  % Add internal psi-functions in W
        % Add coordinates of \bpsi_{0,N_0},...,\bpsi_{0,N'-1}^b in phi_1^b. These come at columns (K-N+2N0+1):(K-N+2(N'-1)+1)
        insertpsi = Gsegment(g1, suppg1, 0:(KminN + 2*(Nprime-1) + 1 + suppg1(end)), KminN + 2*(N0:(Nprime-1)) + 1);
        [W, insertpsi] = expand_cols_smallest(W, insertpsi);
        W( :, KminN + 2*(N0:(Nprime-1)) + 2) = insertpsi;
    end
    if N0 > Nprime % Add internal phi-functions in W
        % Add coordinates of \bphi_{0,N'},...,\bphi_{0,N0-1}^b in phi_1^b. These come at columns (K-N+2N'):(K-N+2(N0-1))
        insertphi = Gsegment(g0, suppg0, 0:(KminN + 2*N0-2 + suppg0(end)), KminN + 2*(Nprime:(N0-1)));
        [W, insertphi] = expand_cols_smallest(W, insertphi);
        W( :, KminN + 2*(Nprime:(N0-1)) + 1) = insertphi;
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