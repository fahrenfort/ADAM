% demo code for MCCA

disp('MCCA applied to real EEG data - no preprocessing');

clear
if exist('../../DATA/MCCA_EXAMPLE_DATA/EEGVolume.mat', 'file')
    load ('../../DATA/MCCA_EXAMPLE_DATA/EEGVolume.mat'); % loads X badchannels eogchannels fs
else
    error('Download data file https://www.parralab.org/isc/EEGVolume.mat');
end

sr=fs;

% Data are organized as a 3D matrix with dimensions time, channel, set
% (subject). These seem to be raw EEG, time-aligned across subjects but
% otherwise unprocessed. In particular each channel has a large offset and
% drift.  We remove the offsets but not the drifts.

[nsamples,nchans,nsets]=size(X);
for iSet=1:nsets
    X(:,:,iSet)=nt_demean(X(:,:,iSet));
end

x=X(:,:);
figure(1); clf
t=(0:nsamples-1)'/sr;
plot(t,x);
xlabel('time'); 
title('raw EEG (all subjects)');


% MCCA
C=x'*x;
[A,score,AA]=nt_mcca(C,nchans);
z=x*A;

figure(2); clf;
subplot 121
plot(score, '.-'); xlabel('SC'); ylabel('variance');
for k=1:4
    subplot(4,2,k*2);
    plot(z(:,k)); title(['SC ',num2str(k)]); axis tight
    set(gca,'ytick',[]);
end
xlabel('samples');

figure(3); clf;
DSR=40; nt_spect_plot2(nt_resample(z(:,1:200),1,DSR),1024,[],[],sr/DSR)
ylabel('SC');
xlabel('frequency (Hz)');
title('power spectra of first SCs');


% -------------------------------------------------------------------------
function X = preprocess(X,eogchannels,badchannels,fs)
% All the usual EEG preprocessing, except epoching and epoch rejection as
% there are not events to epoch for natural stimuli. duh! Instead, bad data
% is set = 0 in the continuous stream, which makes sense when computing
% covariance matrices but maybe not for other purposes. Bad channels are
% removed for all the indice given in badchannels cell array (each subject
% has its own vector of indice). None are removed if this is set to []. If
% it is set to -1, channels are removed based on outlies in power.

debug = 1;     % turn this on to show data before/after preprocessing. 
kIQD=4;        % multiple of interquartile differences to mark as outliers samples
kIQDp=3;       % multiple of interquartile differences to mark as outliers channels
HPcutoff =0.5; % HP filter cut-off frequequency in Hz

% pick your preferred high-pass filter
[z,p,k]=butter(5,HPcutoff/fs*2,'high'); sos = zp2sos(z,p,k);

[T,D,N]=size(X); 

% if it is not EOG, then it must be EEG channel
eegchannels = setdiff(1:D,eogchannels);

% Preprocess data for all N subjects
for i=1:N
    
    data = X(:,:,i);

    % remove starting offset to avoid filter trancient
    data = data-repmat(data(1,:),T,1);
    
    % show the original data
    if debug, subplot(2,1,1); imagesc((1:T)/fs,1:D,data'); title(['Subject ' num2str(i)]); end

    % high-pass filter
    data = sosfilt(sos,data);          
    
    % regress out eye-movements;
    data = data - data(:,eogchannels) * (data(:,eogchannels)\data);     

    % detect outliers above stdThresh per channel; 
    data(abs(data)>kIQD*repmat(diff(prctile(data,[25 75])),[T 1])) = NaN;
    
    % remove 40ms before and after;
    h=[1; zeros(round(0.04*fs)-1,1)];    
    data = filter(h,1,flipud(filter(h,1,flipud(data))));
    
    % Mark outliers as 0, to avoid NaN coding and to discount noisy channels
    data(isnan(data))=0;

    % Find bad channels based on power ourliers, if not specified "by hand"
    if badchannels{i} == -1, 
        logpower = log(std(data)); Q=prctile(log(std(data(:,eegchannels))),[25 50 75]);
        badchannels{i} = find(logpower-Q(2)>kIQDp*(Q(3)-Q(1)));  
    end
    
    % zero out bad channels
    data(:,badchannels{i})=0; 
    
    % show the result of all this
    if debug, subplot(2,1,2); imagesc((1:T)/fs,1:D,data'); caxis([-100 100]); xlabel('Time (s)'); drawnow; end

    X(:,:,i) = data;
    
end

% remove the eog channels as we have used them already
X = X(:,eegchannels,:);
end
