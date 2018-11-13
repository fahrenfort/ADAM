% Script to reproduce Fig. 4 of:
% de Cheveigné, A., Arzounian, D., Robust detrending, rereferencing,
% outlier detection, and inpainting for multichannel data,(submitted to
% NeuroImage).
%
% This script requires the NoiseTools toolbox
% (http://audition.ens.fr/adc/NoiseTools/). Tested with version
% 29-Nov-2017.

clear
datadir='../../../DATA/figures_deCheveigne_Arzounian_2017/';

%% Load EEG data from file
if 0 % Initial creation of MAT data file from raw source file   
    % Ensure that biosig works:
    if 2 ~=exist('sopen')
        error('download BioSig from http://biosig.sourceforge.net');
    end
    x = fileparts( which('sopen') );
    rmpath(x);
    addpath(x,'-begin');
    % Get data from raw source file:
    fname='/data/meg/arzounian/ADC_DA_140521_p20/ADC_DA_140521_p20_01_calib';
    [p,x]=nt_read_data(fname);
    % Reshape to time x channel:
    x=x'; 
    % Remove channel means
    x=nt_demean(x);
    % Convert from V to mV
    x=x/1000; 
    % Select channel of interest:
    xx=x(:,10);
    % Save relevant data
    sr=p.sr;
    save ([datadir, 'fig_4_data'], 'sr', 'xx', 'fname');
else % Load data from MAT file
    if 2~=exist([datadir,'fig_4_data.mat'],'file')
        error('Download from http://audition.ens.fr/adc/NoiseTools/DATA/figures_deCheveigne_Arzounian_2017/');
    end
    load ([datadir, 'fig_4_data']);
end

%% Superimpose repeated, artificial pulse on real EEG

% Remove last portion of EEG, artifacted:
xx=-xx(1:2500000);
t=(0:size(xx,1)-1)'/sr; % time vector

% Create artificial pulse
nsamples=3*sr; % epoch length
signal=zeros(nsamples,1); % signal baseline
n=sr/2; % pulse length
signal(sr/2+(1:n))=sin(pi*(1:n)'/n); % half sinusoid
tt=(0:nsamples-1)'/sr; % epoch time vector

% Add to signal in repeated epochs
ntrials=200;
HOP=10000; % number of samples between consecutive pulse onsets
for iTrial=1:ntrials
    xx(iTrial*HOP+(1:nsamples),:)=xx(iTrial*HOP+(1:nsamples),:)+0.02*signal;
end 

% Plot overall signal time-course
figure(1); clf; set(gcf,'color', [1 1 1], 'position',[440   585   600   150])
plot(t,xx,'k','linewidth',1);
set(gca,'fontsize',14, 'ytick', [-2 0 2],'xtick',0:200:1200); xlabel('time (s)'); 
xlim([t(1) t(end)]);  ylim([-3 4])
ylabel('amplitude (mV)');
title('raw')

% Save panel
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_4_top');

%% Trial-average

% Cut signal into epochs
xxx=zeros(nsamples,size(xx,2),ntrials);
for iTrial=1:ntrials
    xxx(:,:,iTrial)= xx(iTrial*HOP+(1:nsamples),:);
end
% Subtract mean of each trial
xxx=nt_demean2(xxx);

% Plot the trend of average over epochs
figure(2);clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 350])
subplot 221
% Plot normal linear fit of the average over trials
plot(tt,mean(xxx,3)-nt_detrend(mean(xxx,3),1,[],[],100),':r','linewidth',1); hold on
% Plot robust linear fit of average
plot(tt,mean(xxx,3)-nt_detrend(mean(xxx,3),1,[],[],2),'r','linewidth',2);
% Plot average
plot(tt,mean(xxx,3),'k','linewidth',1);
% Add zero line
plot(tt,tt*0,':k');
% Set axes properties and legend
set(gca,'fontsize',12, 'ytick', .02*[-1 0 1], 'xtick',[0 1 2]); xlabel('time (s)'); ylabel('amplitude (mV)'); xlim([tt(1) tt(end)]);
ylim(.02*[-1 1.2]); 
h=legend('fit','robust fit'); legend boxoff
title('average over trials');

%% Detrending

% Plot the signal detrended by subtracting 'non-robust' fit
subplot 223
plot(tt,nt_detrend(mean(xxx,3),1,[],[],100),'k','linewidth',1);
hold on; plot(tt,tt*0,':k');
set(gca,'fontsize',12, 'ytick', .02*[-1 0 1], 'xtick',[0 1 2]); xlabel('time (s)'); ylabel('amplitude (mV)'); xlim([tt(1) tt(end)]);
ylim(.02*[-1 1.2])
title('detrend');

% Plot the signal detrended by subtracting robust fit
subplot 224
plot(tt,nt_detrend(mean(xxx,3),1,[],[],2),'k','linewidth',1);
hold on; plot(tt,tt*0,':k');
mask=ones(size(tt));
mask(sr/2+(1:n))=nan;
plot(tt,-0.4*mask,'color',0.6*[1 1 1],'linewidth', 4);
set(gca,'fontsize',12, 'ytick', [], 'xtick',[0 1 2]); xlabel('time (s)'); xlim([tt(1) tt(end)]);
ylim(.02*[-1 1.2])
title('robust detrend')

%% Trial-average of high-pass data

% High-pass the continuous signal above 0.3 Hz with 4th-order Butterworth filter
[B,A]=butter(4,0.3/(sr/2),'high'); % Designs the filter
xx=filter(B,A,xx); % Causal filtering

% Re-build 3-D matrix with high-pass filtered signal
xxx=zeros(nsamples,size(xx,2),ntrials);
for iTrial=1:ntrials
    xxx(:,:,iTrial)= xx(iTrial*HOP+(1:nsamples),:);
end

% Plot average of high-pass filtered signal over trials
subplot 222
plot(tt,mean(xxx,3),'k','linewidth',1);
hold on; plot(t,t*0,':k', 'linewidth',0.5);
set(gca,'fontsize',12, 'ytick',[], 'xtick',[0 1 2], 'box','on'); xlabel('time (s)'); xlim([tt(1) tt(end)]);
ylim(.02*[-1 1.2])
title('high-pass 0.3 Hz order=4');

% Save 4 panels
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_4_bottom');

