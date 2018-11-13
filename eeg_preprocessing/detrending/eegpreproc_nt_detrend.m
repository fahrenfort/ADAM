function eegpreproc_nt_detrend(subject_path,epoch_pad,polynom)
% subject-specific function that can be run in parallel on Lisa to preprocess multiple subjects at
% once
%
% requires path to raw data file as first input
%
% this preprocessing function runs robust detrending instead of high-pass filtering; this requires
% one input parameter to specify the epoch length used for detrending; this value is padded on each
% side of the original epoch (e.g. epoch_pad of '100' creates 200+ sec epochs); the other input
% parameter is the polynomial order
%
% parameters are given in strings and internally converted to double, such that it works on lisa
% after compiling
%
% other preprocessing steps are a bit idiosyncratic for this study: bad channels to interpolate were
% already selected in previous analyses and saved into matfiles that are located in the input
% folder of the raw data; these are loaded and interpolated; error trials and missings are removed;
% muscle artifacts are not removed, this is done with ADAM; ICA is not performed.

%
cd(subject_path)
fprintf('Preprocess subject %s\n',subject_path)

% the data were already re-referenced to average earlobes
fn = dir('*reref.mat');
load bad_chans

epoch_pad = str2double(epoch_pad);
polynom = str2double(polynom);

% preproc ingredients

% let's pick all short trials, also the short trials from the mixed block that contained long trials
% as well; this difference is not very important for the illustration of filtering/detrending, and
% increases trial count, which is preferred in this case

triggers={
    'face_short'         { '11' };
    'house_short'        { '12' };
    'letter_short'       { '13' };
%     'face_long'          { '21' };
%     'house_long'         { '22' };
%     'letter_long'        { '23' };
    'face_shortmixed'         { '31' };
    'house_shortmixed'        { '32' };
    'letter_shortmixed'       { '33' };
    };

cfg = [];
cfg.writdir         = pwd;
cfg.triggers        = triggers;
cfg.nchan           = 64;
% two separate time windows: one for the actual epoch of interest, one for detrending which is this
% same epoch plus some time padded to both sides, specified by epoch_pad; this parameter controls
% the goodness of fit of the polynomial order; currently we pick 100-200 range according to examples
% given in de Cheveign? & Arzounian 2018.
cfg.epochtime       = [-1.5 6]; % 250 ms cue, 3000 ms retention, RT of ~1500 ms; 1.5 sec on each side gives this epoch
cfg.detrendtime     = cfg.epochtime + [-1*epoch_pad epoch_pad];
cfg.polynom         = polynom;
cfg.badchans        = bad_chansidx;
cfg.filename        = fn(1).name;
cfg.outfilename     = [fn(1).name(1:4) '_temporientfollowup_detrend' num2str(polynom) '_epoch' num2str(epoch_pad) '.set' ];

% PREPROCESSING WITH ROBUST DETRENDING
% unpack cfg
v2struct(cfg);

% load data

fprintf('Loading data...\n');
load(filename)
% original sampling rate was 2048; in the re-referencing stage resampled to 512; now further
% downsampled to 256 (and during decoding to 32 Hz)
EEG = pop_resample(EEG,256);

% remove EOG for good
EEG = pop_select(EEG,'nochannel',[nchan+1 nchan+2]);

% mirror-pad edges

EEG.data = padarray(EEG.data,[0 epoch_pad*EEG.srate],'both','symmetric');
tdiff = EEG.times(2)-EEG.times(1);
EEG.times(end+1:end+2*epoch_pad*EEG.srate) = EEG.times(end)+tdiff:tdiff:EEG.times(end)+tdiff*(2*epoch_pad*EEG.srate);
EEG.pnts = length(EEG.times);
EEG.xmin = EEG.times(1);
EEG.xmax = EEG.times(end);

for evi=1:length(EEG.event)
    EEG.event(evi).latency = EEG.event(evi).latency+EEG.srate*epoch_pad;
    EEG.urevent(evi).latency = EEG.urevent(evi).latency+EEG.srate*epoch_pad;
end

% Epoch

EEG = pop_epoch( EEG, [triggers{:,2}],detrendtime);

% find error trials and missings

errortrials=zeros(1,EEG.trials);
missing = errortrials;

for ei=1:EEG.trials
    try
        
        EEG.epoch(ei).trialnum = ei;
        latencies = cell2mat(EEG.epoch(ei).eventlatency);
        events = cell2mat(EEG.epoch(ei).eventtype);
        
        search = ismember(events,50:56); search=find(search); search(latencies(search)<0)=[]; 
        cues = ismember(events,11:43); cues=find(cues); cues(latencies(cues)<0)=[]; 
        response = ismember(events,[512 1024]); response=find(response); response(latencies(response)<0)=[]; 

        EEG.epoch(ei).rt = latencies(response(1))-latencies(cues(1));
        if length(cues)>1
            EEG.epoch(ei).next = latencies(cues(2))-latencies(cues(1));
        else
            EEG.epoch(ei).next = 9999;
        end
        if EEG.epoch(ei).rt > 3000+6267
            missing(ei)=1;
        end
        
        % task was an absent/present judgment and below triggers correspond to this mapping
        if search(1)==50 && response(1)==512
            errortrials(ei)=1;
        elseif search(1)>50 && response(1)==1024
            errortrials(ei)=1;
        end
        
    catch me
        missing(ei)=1;
        errmess{ei}=me;
    end
end

% remove errors missings

EEG = pop_select(EEG,'notrial',find(errortrials+missing));

% use NoiseTools for robust detrending
% toolbox can downloaded from: audition.ens.fr/adc/NoiseTools

tic
fprintf('Robust detrending @ %ith order polynomial\n',polynom);

filt_dat = zeros(size(EEG.data));
for ei=1:EEG.trials
    
    wt = ones(EEG.pnts,1);
    
    % mask is pre-set to exclude [0 - (RT+750ms)] interval from the fit
    tidx = dsearchn(EEG.times',[0 EEG.epoch(ei).rt+750]')';
    wt(tidx(1):tidx(2))=0;
    
    [tmp,w1] = nt_detrend(EEG.data(:,:,ei)',1,repmat(wt,[1 length(EEG.chanlocs)])); % start with 1st order
    if polynom>1
        [x2, ~] = nt_detrend(tmp,polynom,w1); % then nth order with mask of previous step
        filt_dat(:,:,ei) = x2';
    else
        filt_dat(:,:,ei) = tmp';
    end
end
toc

%
EEG.data = filt_dat;
clear filt_dat

% narrow epoch

EEG = pop_select( EEG, 'time',epochtime);

% interpolate bad channels

if sum(bad_chansidx)>0
    EEG = eeg_interp(EEG,bad_chansidx);
end


pop_saveset(EEG,'filename',outfilename);

