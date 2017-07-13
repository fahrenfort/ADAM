%------------------------ toolboxes ------------------------%

% path definitions
ft_path = fullfile('C:','matlab_toolboxes','fieldtrip-20170704'); % Mac starts with file separator, example: fullfile(filesep,'Users','Johannes','matlab_toolboxes','fieldtrip-20170704');
eeglab_path = fullfile('C:','matlab_toolboxes','eeglab14_1_1b');
adam_path = fullfile('C:','matlab_toolboxes','ADAM-master');

% FT
if (exist(ft_path,'dir') == 7) && (~isdeployed)
    addpath(ft_path,'-begin');
    ft_defaults; % find default options and files
    disp('FIELDTRIP IS ALIVE');
elseif ~isdeployed
    disp(['WARNING, CANNOT FIND FIELDTRIP TOOLBOX AT ' ft_path ', CHECK PATHS IN startup.m']);
end

% EEGLAB
if (exist(eeglab_path,'dir') == 7) && (~isdeployed)
    curdir = pwd;
    cd(eeglab_path);
    eeglab;
    cd(curdir);
    disp('EEGLAB IS ALIVE');
elseif ~isdeployed
    disp(['WARNING, CANNOT FIND EEGLAB TOOLBOX AT ' eeglab_path ', CHECK PATHS IN startup.m']);
end
if ~isdeployed
    % remove conflicting paths
    rmpath(genpath(fullfile(eeglab_path,'external','fieldtrip-partial')));
    % create eeglabexefolder function to replace eeglab's internal
    % function and put options files in eeglab root
    tmpf = which('eeglabexefolder.m');
    if ~isempty(tmpf)
        if isempty(which('eeglabexefolder_original.m'))
            movefile(tmpf,[ tmpf(1:end-2) '_original.m']);
        end
        fid = fopen(tmpf,'w');
        fprintf(fid, 'function eeglabdir = eeglabexefolder\n');
        fprintf(fid, '%% This is to replace eeglabs native function, so it knows where to find the option files after compiling.\n');
        fprintf(fid, 'eeglabdir = ''%s'';\n', eeglab_path);
        fclose(fid);
        tmpf = which('eeg_optionsbackup.m');
        copyfile(tmpf, fullfile(eeglabexefolder, 'eeg_optionsbackup.txt'));
        tmpf = which('eeg_options.m');
        copyfile(tmpf, fullfile(eeglabexefolder, 'eeg_options.txt')); 
        disp(['EEGLAB: exefolder is here ' eeglabexefolder]);
    end
end

% ADAM decoding toolbox
if (exist(adam_path,'dir') == 7) && (~isdeployed)
    addpath(genpath(adam_path),'-begin');
    % create findcapfile function to return location of capfile
    tmpf = which('standard-10-5-cap385.elp');
    if ~isempty(tmpf)
        fid = fopen(fullfile(fileparts(tmpf),'findcapfile.m'),'w');
        fprintf(fid, 'function capfile = findcapfile\n');
        fprintf(fid, '%% This is to be able to find the capfile for elecrode lookup after compiling.\n');
        fprintf(fid, 'capfile = ''%s'';\n', tmpf);
        fclose(fid);
    end
    disp('ADAM IS ALIVE');
    disp(['ADAM:   capfile is here ' findcapfile]);
else
    disp(['WARNING, CANNOT FIND ADAM TOOLBOX AT ' adam_path ', CHECK PATHS IN startup.m']);
end

% clear rubish
close;
clear;