




function filter=getDBfilter(vm, type)
    
    %cur_dir =  fileparts(mfilename('fullpath'));
    %
    %if (type == 0)
    %    filename = sprintf('%s/var/DB%d.mat', cur_dir, vm);
    %else
    %    filename = sprintf('%s/var/sym%d.mat', cur_dir, vm);
    %end
    %
    %if (exist(filename) == 2);
    %    load(filename);
    %else
        filter = liftingfactortho(vm, type, 1);
        
    %    if (exist(strcat(cur_dir, '/var')) ~= 7)
    %        mkdir(strcat(cur_dir,'/var'));
    %    end
    %    
    %    save(filename, 'filter');
    %end
end





