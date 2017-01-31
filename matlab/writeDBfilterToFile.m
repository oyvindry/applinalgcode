function writeDBfilterToFile(vm_max)
    
    dest = 'var';
    
    if (exist(dest) ~= 7)
        mkdir(dest);
    end
    
    for vm = 2:vm_max
        filter = liftingfactortho(vm, 0);
        save(sprintf('%s/DB%d.mat', dest, vm), 'filter');
        filter = liftingfactortho(vm, 1);
        save(sprintf('%s/sym%d.mat', dest, vm), 'filter');
    end
end
