function filter = getDBfilter(vm, type)
    
    dest = 'var';
    if (type == 0)
        filename = sprintf('%s/DB%d.mat', dest, vm);
    else
        filename = sprintf('%s/sym%d.mat', dest, vm);
    end

    if (exist(filename) == 2);
        load(filename);
    else
        filter = liftingfactortho(vm, type);
        
        if (exist(dest) ~= 7)
            mkdir(dest);
        end
        
        save(filename, 'filter');
    end
end
