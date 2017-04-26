%------------ toolboxes ------------------------------%

% path definitions
ft_path = sprintf('/home/johannes/matlab_toolboxes/fieldtrip-20150318');
%ft_path = sprintf('/home/johannes/matlab_toolboxes/fieldtrip-20161122');
eeglab_path = sprintf('/home/johannes/matlab_toolboxes/eeglab13_2_2b');
%eeglab_path = sprintf('/home/johannes/matlab_toolboxes/eeglab13_6_5b');
adam_path = sprintf('/home/johannes/Dropbox/matlab_scripts');

% FT
if (exist(ft_path,'dir') == 7) && (~isdeployed)
    addpath(ft_path,'-begin');
end

% EEGLAB
if (exist(eeglab_path,'dir') == 7) && (~isdeployed)
     addpath(eeglab_path);
     eeglab;
     
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
     end
     tmpf = which('eeg_optionsbackup.m');
     copyfile(tmpf, fullfile(eeglabexefolder, 'eeg_optionsbackup.txt'));
     tmpf = which('eeg_options.m');
     copyfile(tmpf, fullfile(eeglabexefolder, 'eeg_options.txt'));
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
end

% find default options and files
ft_defaults;
disp(['EEGLAB: exefolder is here ' eeglabexefolder]);
disp(['ADAM:   capfile is here ' findcapfile]);

% clear rubish
clear;