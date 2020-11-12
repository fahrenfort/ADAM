%% this code is meant to obtain trigger codes from somewhere else and save them
logdatadir = '/Users/VU-MBP/Dropbox/Work/Experimenten/WM_debunk_EEG/LOGFILES';
logfiles = {'S001_WM_debunk_2015_Jun_18_0953', 'S001_WM_debunk_2015_Jun_23_1007'}; 
% logstruct = readPsychoPy(logdatadir,logfiles); readPsychoPy doesn't work so well, use proprietary
nFiles = numel(logfiles);
for cFile = 1:nFiles
    % read data
    triggerstruct.(logfiles{cFile}) = readtable([logdatadir filesep logfiles{cFile} '.csv']);
end
save('/Users/VU-MBP/Dropbox/Work/Experimenten/WM_debunk_EEG/info/newtriggerinfo1.mat', 'triggerstruct');
% you can load newtriggerinfo1 later to import triggers or something
