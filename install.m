% This install file tries to detect your 'startup.m' file and add the requiered
% wavlib paths as default. It also add these paths to the current session.
% If this install file fails, you will have to find the 'startup.m' file 
% and add the requiered paths manually. 

% IMPROVE
% In the next version of this file, I should consider forking out the part of
% the code which  

dependent_paths = {'matlab', 'sounds', 'images'};

us_root = userpath();
startup_file = strcat(us_root(1:end-1), '/startup.m'); 

if (exist(us_root(1:end-1)) == 7)
    
    fID = fopen(startup_file, 'a');

    fprintf(fID, strcat(...
                 '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', ...
                 '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', ...
                 '%% These paths have automatically been added by the\n', ... 
                 '%% install.m belonging to wavlib\n'));
    
    for wavlib_path = dependent_paths 
        fprintf(fID, 'addpath %s/%s\n', pwd(), wavlib_path{1});
    end
    fprintf(fID, strcat( ... 
                 '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', ...
                 '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'));
    fclose(fID);
    fprintf('wavlib have added the following paths:\n');
    for wavlib_path = dependent_paths 
        fprintf('addpath %s/%s\n', pwd(), wavlib_path{1});
        addpath(sprintf('%s/%s', pwd(), wavlib_path{1}))
    end
    fprintf(strcat(...
                 'to the file ', ... 
                 ' %s\n'), startup_file);
else 
    fprintf('wavlib was unable to find the file startup.m\n');        
    fprintf('in your `userpath()`.');
    fprintf('Adding the following paths to the current session');
    for wavlib_path = dependent_paths 
        fprintf('addpath %s/%s\n', pwd(), wavlib_path{1});
        addpath(sprintf('%s/%s', pwd(), wavlib_path{1}));
    end
end

