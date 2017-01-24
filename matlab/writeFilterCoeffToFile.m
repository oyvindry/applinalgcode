function writeFilterCoeffToFile( W, filename );
    
    N = size(W,2)/2;
    
    S = N+1;
    
    fid = fopen(filename, 'w');
    for i=1:N
        for k=1:(S+2*(i-1))
            
            fprintf(fid, '%e %e\n', W(k, 2*i-1), W(k, 2*i) );
            fprintf(fid, '\n');
        end
        fprintf(fid, '\n\n');
    end
    
    fclose(fid);
    

