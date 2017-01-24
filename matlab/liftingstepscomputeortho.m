function filters = liftingstepscomputeortho(h0, h1)
    
    %global wavlib_lambdas wavlib_alpha wavlib_beta h0 h1
    
    stepnr=1;
    len1=length(h0)/2; len2=len1;
    filters.lambdas=zeros(len1+1,2);
    if mod(len1,2)==0
        h00=h0(1:2:length(h0));
        h01=h0(2:2:length(h0));
        h10=h1(1:2:length(h1));
        h11=h1(2:2:length(h1));
  
        lambda1=-h00(1)/h10(1);
        h00=h00+lambda1*h10; 
        h01=h01+lambda1*h11;
        start1=2; end1=len1; len1=len1-1; start2=1; end2=len2;
        filters.lambdas(stepnr,:)=[lambda1 0];
    else
        h00=h0(2:2:length(h0));
        h01=h0(1:2:length(h0));
        h10=h1(2:2:length(h1));
        h11=h1(1:2:length(h1));
    
        lambda1=-h10(len1)/h00(len1); 
        h10=h10+lambda1*h00; 
        h11=h11+lambda1*h01;
        start2=1; end2=len2-1; len2=len2-1; start1=1; end1=len1;
        filters.lambdas(stepnr,:)=[0 lambda1];
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
        filters.lambdas(stepnr,:)=[lambda1 lambda2];
        stepnr=stepnr+1;
    
        %[h00 h01; h10 h11]
        %conv(h00,h11)-conv(h10,h01)
    end
  
    % Add the final lifting, and alpha,beta
    filters.alpha=sum(h00);
    filters.beta=sum(h11);
    lastlift=-sum(h01)/filters.beta;
    if mod(length(h0)/2,2)==0
        filters.lambdas(stepnr,:)=[0 lastlift];
    else
        filters.lambdas(stepnr,:)=[lastlift 0];
    end
    %[h00 h01; h10 h11]
