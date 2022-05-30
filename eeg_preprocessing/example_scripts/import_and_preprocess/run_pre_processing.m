%% trigger labels and epoch settings
clear;

% file names
filenames= {'PP001.set'     'PP002.set' ...
            'PP003.set'     'PP004.set' ...
            'PP005.set'     'PP006.set' };

% all cues
cue_face = [1:3 11:13 21:23 31:33];
cue_house = [4:6 14:16 24:26 34:36];
cue_letter = [7:9 17:19 27:29 37:39];

% all stimuli
stim_face = [1 4 7 11 14 17 21 24 27 31 34 37];
stim_house = [2 5 8 12 15 18 22 25 28 32 35 38];
stim_letter = [3 6 9 13 16 19 23 26 29 33 36 39];

% targets and non-targets in prediction session
upright = 1:9;                        % non-target (do not press)
upsidedown = 11:19;                   % target (press space bar)
predicted = [1 5 9 11 16 19];         % targets and non-targets (no relation to button press)
nonpredicted = [2:4 6:8 12:14 16:18]; % targets and non-targets (no relation to button press)
triggers = [predicted nonpredicted];

% epoch time period
eeglab_epochtime = [-2.5 1];

%% pipeline to follow:
% - import into EEGLAB
% - highpass filter of .01 Hz
% - epoch (-2.5 1) ms
% - remove blinks using ICA
%
% During analysis:
% - reject trials with muscle artifacts in a window from 0-500 ms (large deviating power in the 110-140 Hz range)
% - baseline correction on (-2300,-2000) ms -> the fixation period 
% - lots of decoding on lots of stuff

%% read data into eeg lab
files2do = {};
basedir = '/Volumes/backup/Predictive_EEG/RAW_logfiles_bdfs/BDFs';
save2dir = '/Volumes/backup/Predictive_EEG/EEGLAB_DATA/unfiltered_unepoched';
filz=dir([basedir filesep '*.bdf']); % find all set files for each subject
for cFiles=29:numel(filz)
    filename = [basedir filesep filz(cFiles).name];
    if sum(strncmp(filz(cFiles).name,files2do,numel(filz(cFiles).name)-4))>0 || isempty(files2do)
        EEG = readbdf_into_eeglab_vu(filename);
        % downsample a bit, 512 is really not necessary
        EEG = pop_resample(EEG, 256);
        % get the name
        [~, fname, ~ ] = fileparts(filename);
        EEG = pop_saveset(EEG, 'filename',strrep(fname,'.','_'),'filepath',save2dir);
        events = [EEG.event(:).type];
        trig_events = sum(ismember(events, triggers));
    end
end

%% filter and epoch data on server
sourcedir =  '$HOME/predictive/EEGLAB_DATA/unfiltered_unepoched';
targetdir =  '$HOME/predictive/EEGLAB_DATA/highpass_.01';

triglabels = prediction;
filter_eeg_and_epoch(sourcedir,filenames,targetdir,0.01,0,eeglab_epochtime(1),eeglab_epochtime(2),cond_string(triggers));

%% run ICA
filepath = '$HOME/predictive/EEGLAB_DATA/highpass_.01';
outpath = '$HOME/predictive/EEGLAB_DATA/IC';
for cFiles=1:numel(filenames)
    filename = filenames{cFiles};
    compute_ICs_new(filepath,filename,outpath,'no','no','no','yes');
end

%% Remove ICs 
sourcedir = '$HOME/predictive/EEGLAB_DATA/ICA';
targetdir = '$HOME/predictive/EEGLAB_DATA/noblinks';
for cFiles=1:numel(filenames)
    remove_ICs(sourcedir,filenames,targetdir,'0,0,1,0');
end

% %% Remove eyeblinks RUN IN MATLAB ON LISA
% clear;
% sourcedir = '/home/johannes/N2pc_Eimer2/EEGLAB_DATA/noblinks';
% targetdir = '/home/johannes/N2pc_Eimer2/EEGLAB_DATA/noblinks_noeye';
% filz=dir([sourcedir filesep '*.set']); % find all set files
% for cFiles=1:numel(filz)
%     filename = filz(cFiles).name;
%     fullfile = [sourcedir filesep filename];
%     EEG = pop_loadset('filename',filename,'filepath',sourcedir);
%     events = pop_exportepoch(EEG);
%     % make sure you get the labels for the ocular channels right
%     channelnames = {EEG.chanlocs(:).labels};
%     timewindow = [0 500]; % time window to look at
%     % compute pure horizontal eye movement channel
%     EEG2 = EEG;
%     EEG2.data(ismember(channelnames,{'HEOG_R'}),:,:) = EEG2.data(ismember(channelnames,{'HEOG_R'}),:,:) - EEG2.data(ismember(channelnames,{'HEOG_L'}),:,:);
%     % find eye-movements
%     rejected_trials = pop_artstepX(EEG2,'channel',find(ismember(channelnames,{'HEOG_R'})),'flag',1,'threshold',30,'Twindow', timewindow, 'Windowsize',  100, 'Windowstep',  50 );
%     clear EEG2;
%     rejected_trials = logical(rejected_trials);
%     EEG.reject.rejmanual = rejected_trials;
%     % add 10000
%     events(rejected_trials) = events(rejected_trials) + 10000;
%     EEG = pop_importepoch(EEG,events,{'eventtype'},'typefield','eventtype');
%     % save
%     pop_saveset(EEG, 'filename',filename,'filepath',targetdir);
%     filenames{cFiles,1} = filename;
%     percentage_eye_movements(cFiles,1) = sum(rejected_trials)/numel(rejected_trials);
% end

% %% EXAMPLE Recode incorrect responses by adding 1000 for missed and incorrect trials
% sourcedir = '/Users/VU-MBP/Dropbox/Work/Experimenten/N2pc_Eimer2/EEGLAB_DATA/highpass_.1';
% targetdir = '/Users/VU-MBP/Dropbox/Work/Experimenten/N2pc_Eimer2/EEGLAB_DATA/recoded';
% filz = dir([sourcedir filesep '*.set']);
% clear percentage;
% for cFiles = 1:numel(filz)
%     fname = filz(cFiles).name;
%     disp(fname);
%     EEG = pop_loadset('filename',fname,'filepath',sourcedir);
%     n = 0;
%     for cEpoch= 1 : numel(EEG.epoch)
%         % then recode incorrect response trials by adding 1000
%         if sum(strcmp(EEG.epoch(cEpoch).eventtype,'100')) == 0 % 100 are the correct responses, so this identifies incorrect or absent responses
%             %disp(['recoding event ' num2str(cEpoch)]);
%             indxoftrig = cell2mat(EEG.epoch(cEpoch).eventlatency) == 0;
%             if ~isempty(indxoftrig)
%                 n = n + 1;
%                 eventcode = string2double(EEG.event(EEG.epoch(cEpoch).event(indxoftrig)).type);
%                 EEG.event(EEG.epoch(cEpoch).event(indxoftrig)).type = num2str(eventcode+1000);
%             end
%         end
%     end
%     disp(['recoded ' num2str(n) ' events.']);
%     percentage(cFiles) = n/numel(EEG.epoch);
%     EEG = eeg_checkset(EEG,'eventconsistency');
%     EEG = pop_saveset(EEG, 'filename',fname,'filepath',targetdir);
% end

% %% EXAMPLE count events for check
% for cEpoch = 1 : numel(EEG.epoch)
%     indxoftrig = cell2mat(EEG.epoch(cEpoch).eventlatency) == 0;
%     event_val(cEpoch) = string2double(EEG.event(EEG.epoch(cEpoch).event(indxoftrig)).type);
% end