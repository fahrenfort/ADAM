%% this code is meant to get the correct trigger codes
logdatadir = '/Users/VU-MBP/Dropbox/Work/Experimenten/WM_debunk_EEG/RAW_logfiles_bdfs/logfiles';
logfiles = {'logfile1', 'logfile2'}; 
% read data
logstruct = readPsychoPy(logdatadir,logfiles);
nSes = numel(logfiles);
% fix trigger info
for cSes = 1:nSes
    data = logstruct.(logfiles{cSes});
    % here some computations to compute the actual conditions
    % the result you put in actual_triggers
    actual_triggers = pop_exportepoch(EEG); % pop_exportepoch obtains the event values of epoched data
    actual_triggers = actual_triggers + 6; % do something to triggers, here add 6 for illustration purposes
    triggerstruct.(logfiles{cSes}) = actual_triggers;
end
save('/Users/VU-MBP/Dropbox/Work/Experimenten/WM_debunk_EEG/info/newtriggerinfo.mat', 'triggerstruct');
% copy the above file over to lisa in the sourcedir defined in the first
% line of the code in the cell below

%% this code is meant to update the EEGlab files
% easiest to run this on Sara inside matlab
cd ~;
sourcedir =  'WM_debunk_EEG/EEGLAB_DATA/highpass_.1_epoched';
save2dir = 'WM_debunk_EEG/EEGLAB_DATA/highpass_.1_epoched_fixed_codes';
filz={'file1.set' 'file2.set'};
load([sourcedir filesep 'newtriggerinfo.mat']);
for cFiles=1:numel(filz)
    filename = filz(cFiles).name;
    fullfile = [sourcedir filesep filename];
    EEG = pop_loadset('filename',filename,'filepath',sourcedir);
    newevents = triggerstruct.(['log' filename]);
    EEG = pop_importepoch(EEG,newevents,{'eventtype'},'typefield','eventtype');
    EEG = pop_saveset(EEG,'filename',filename,'filepath',save2dir);
    disp(['saving ' filename]);
end