% Script to reproduce Fig. 8 of:
% de Cheveigné, A., Arzounian, D., Robust detrending, rereferencing,
% outlier detection, and inpainting for multichannel data,(submitted to
% NeuroImage).
%
% This script requires the NoiseTools toolbox
% (http://audition.ens.fr/adc/NoiseTools/). Tested with version
% 29-Nov-2017.

%% Load EEG data from file
clear;
datadir='../../../DATA/figures_deCheveigne_Arzounian_2017/';

if 0 % Initial creation of MAT data file from raw source file   
    % Ensure that biosig works:
    if 2 ~=exist('sopen')
        error('download BioSig from http://biosig.sourceforge.net');
    end
    x = fileparts( which('sopen') );
    rmpath(x);
    addpath(x,'-begin');    
    % Get data from raw source file:
    fname='/data/meg/theoldmanandthesea/eeg/mc/MC_aespa_speech_45.mat';
    [p,x]=nt_read_data(fname); 
    % Reshape to time x channel
    x=x'; 
    % Remove non-EEG channels
    x=x(:,1:128); 
    % Keep only a few seconds of data:
    x=x(0000+(1:6000),:);
    % Convert from V to mV
    x=x/1000; 
    % Sampling rate:
    sr=512; %Hz
    % Time vector
    t=(0:size(x,1)-1)/sr;
    % Save relevanta variables:
    save ([datadir, 'fig_8_data'], 'x', 'sr', 't');
    
else % Load data from MAT file
    if 2~=exist([datadir,'fig_8_data.mat'],'file')
        error('Download from http://audition.ens.fr/adc/NoiseTools/DATA/figures_deCheveigne_Arzounian_2017/');
    end
    load ([datadir,'fig_8_data']);
end

%% Preprocess signals:

% Detrend signals:
x=nt_demean(x); % subtract mean of each channel
[y,w]=nt_detrend(x,3,[],[],3); % get weight vector for robust fit of 3rd order polynomial
[y,w]=nt_detrend(x,10,w,[],3); % subtract 10th order polynomial

% Detect outliers automatically and replace them by signal interpolation:
[ww,yy]=nt_outliers(y,[],1); drawnow;

%% Plot results:

% Prepare figure:
figure(1); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 600])

% Plot estimated weight matrix at the top:
subplot 311; 
nt_imagescc(ww'*2-1); ylabel('channel'); 
% Set axes and colormap properties:
title ('weight (estimated)');
set(gca,'fontsize',14, 'xticklabel',[]);
c=colormap; colormap([repmat(c(1,:),32,1);repmat(c(end,:),32,1)]); 

% Plot raw, detrended and inpainted signals of channel 35 in the middle:
subplot 312;
ch=35; % channel
plot(t,x(:,ch),'k'); hold on; % raw signal
plot(t,y(:,ch)-2 ,'g'); hold on; % detrended signal
plot(t,yy(:,ch)-4,'r'); hold on; % inpainted signal
% Set axes properties and legend:
set(gca,'fontsize',14);
xlim([t(1) t(end)]); ylim([-20 10])
ylabel('amplitude (mV)');
title ('channel #35'); 
set(gca,'fontsize',14, 'xticklabel',[]);
legend('raw', 'detrended', 'fixed','location', 'southeast'); legend boxoff

% Plot raw, detrended and inpainted signals of channel 39 in the middle:
subplot 313;
ch=39; % channel 
plot(t,x(:,ch)+5,'k'); hold on; % raw signal
plot(t,y(:,ch) ,'g'); hold on; % detrended signal
plot(t,yy(:,ch)-5,'r'); hold on; % inpainted signal
% Set axes properties and legend:
set(gca,'fontsize',14);
xlim([t(1) t(end)]); ylim([-20 10])
xlabel('time (s)'); ylabel('amplitude (mV)');
title ('channel #39'); 
legend('raw', 'detrended', 'fixed','location', 'southeast'); legend boxoff

% Save 3 panels:
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_8')
