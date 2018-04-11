function [W,Wtilde,A_pre,A_pre_inv,Atilde_pre,Atilde_pre_inv]=bw_compute_left_biortho(g0,N,h0,Ntilde,debugMode)
    % bw_compute_left_biortho(sym([1/4; 1/2; 1/4]), 2, sym([-1/4; 1/2; 3/2; 1/2; -1/4]), 2);
    % general: 
    %   [h0,h1,g0,g1]=compute_spline_filters(N, Ntilde); % internal function
    %   bw_compute_left_biortho(g0,N,h0,Ntilde);

    Nprime = max(N,Ntilde);
    R = (length(g0)-1)/2; Rtilde = (length(h0)-1)/2;
    L  = -R; Ltilde = -Rtilde;
    
    K = -L; Ktilde = K + Ntilde - N;
    if Ktilde < -Ltilde
        Ktilde = -Ltilde;
        K = Ktilde + N - Ntilde;
    end
    
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
    Gup = Geven(g0,L:R,(-R+1):(K-1),(-2*R+2):2:(2*K-2));
    Gdown = Geven(g0,L:R,K:(2*K+R-2),(-2*R+2):2:(2*K-2));
    
    Gtildeup = Geven(h0,Ltilde:Rtilde,(-Rtilde+1):(Ktilde-1),(-2*Rtilde+2):2:(2*Ktilde-2));
    Gtildedown = Geven(h0,Ltilde:Rtilde,Ktilde:(2*Ktilde+Rtilde-2),(-2*Rtilde+2):2:(2*Ktilde-2));
    
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
    
    % Make phitilde staggered
    [Lmatr,Umatr] = lu((flipud(Ctilde((Rtilde+Ktilde-Ntilde):end,:)))');
    invP = flipud(Lmatr');
    P = inv(invP);
    Xtildee = inv(P)*Xtildee*P;
    Ztildee = Ztildee*P;

    % Address when Ntilde and N are different
    if Ntilde > N    
        XG = Geven(g0,L:R,(2*K+L):(Ktilde-1),(2*K):2:(2*Ktilde-2));
        ZG = Geven(g0,L:R,Ktilde:(2*Ktilde+R-2),(2*K):2:(2*Ktilde-2));
                   
        Xe = [Xe zeros(K+L+N,Ntilde-N); Ze(1:(Ntilde-N),:) XG];
        [Ze,ZG] = expand_cols_smallest(Ze((Ntilde-N+1):end,:), ZG);
        Ze = [Ze ZG];
    elseif N > Ntilde
        XG = Geven(h0,Ltilde:Rtilde,(2*Ktilde+Ltilde):(K-1),(2*Ktilde):2:(2*K-2));
        ZG = Geven(h0,Ltilde:Rtilde,K:(2*K+Rtilde-2),(2*Ktilde):2:(2*K-2));
                
        Xtildee = [Xtildee zeros(Ktilde+Ltilde+Ntilde,N-Ntilde); Ztildee(1:(N-Ntilde),:) XG];
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
    
    % psi-funksjoner

    % Construct Xo
    newcols = Geven(g0,L:R,0:(2*T + R + K - N),(2*Nprime):2:(2*T));
    [Rmatr, newcols] = expand_cols_smallest([Xe;Ze], newcols);
    G = [Rmatr newcols];
    
    newcols = Geven(h0,Ltilde:Rtilde,0:(2*T + Rtilde + K - N),(2*Nprime):2:(2*T));
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
    newcols = Geven(h0,Ltilde:Rtilde,0:(2*Ttilde + Rtilde + K - N),(2*Nprime):2:(2*Ttilde));
    [Rmatr, newcols] = expand_cols_smallest([Xtildee;Ztildee], newcols);
    Gtilde = [Rmatr newcols];
    
    newcols = Geven(g0,L:R,0:(2*Ttilde + R + K - N),(2*Nprime):2:(2*Ttilde));
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

    % testcode
    Xe_test = [Xe;Ze]; Xtildee_test=[Xtildee; Ztildee];
    [Xe_test, Xo_test]=expand_cols_smallest(Xe_test, Xo);
    [Xtildee_test, Xtildeo_test]=expand_cols_smallest(Xtildee_test, Xtildeo);
    W = [Xe_test Xo_test];
    W(:,1:2:end) = Xe_test; W(:,2:2:end) = Xo_test; % TODO: adjust placement
    Wtilde = [Xtildee_test Xtildeo_test];
    Wtilde(:,1:2:end) = Xtildee_test; Wtilde(:,2:2:end) = Xtildeo_test; % TODO: adjust placement
    [W, Wtilde]=expand_cols_smallest(W, Wtilde);
    
    W % only dyadic fractions?
    Wtilde % only dyadic fractions?
    W'*Wtilde % Should be identity matrix
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