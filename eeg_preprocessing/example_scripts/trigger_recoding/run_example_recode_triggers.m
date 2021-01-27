%% this code is meant to import and update EEGLAB trigger codes at time zerocd ~;
sourcedir =  '/Users/VU-MBP/Dropbox/Work/Experimenten/WM_debunk_EEG/EEGLAB_DATA';
save2dir = '/Users/VU-MBP/Dropbox/Work/Experimenten/WM_debunk_EEG/EEGLAB_DATA_NEWTRIGGERS';
filz = dir([sourcedir filesep '*.set']); % {'file1.set' 'file2.set' 'file3.set'};
filenames = {filz(:).name};
for cFiles=1:numel(filz)
    filename = filenames(cFiles);
    fullfile = [sourcedir filesep filename];
    EEG = pop_loadset('filename',filename,'filepath',sourcedir);
    actual_triggers = pop_exportepoch(EEG); % pop_exportepoch obtains the event values of epoched data
    new_triggers = actual_triggers + 6; % do something to triggers, here add 6 for illustration purposes
    EEG = pop_importepoch(EEG,new_triggers,{'eventtype'},'typefield','eventtype');
    EEG = pop_saveset(EEG,'filename',filename,'filepath',save2dir);
    disp(['saving ' filename]);
end