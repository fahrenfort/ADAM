% Script to reproduce Fig. 9 of:
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
    % Remove non-EEG channels:
    x=x(:,1:128); 
    % Keep only a few seconds of data:
    x=x(2000+(1:2000),:);
    % Convert from V to mV:
    x=x/1000;
    % Sampling rate:
    sr=512; %Hz
    % Time vector:
    t=(0:size(x,1)-1)/sr;
    % Save relevant variables:
    save ([datadir,'fig_9_data'], 'x', 'sr', 't');
else % Load data from MAT file
    if 2~=exist([datadir,'fig_9_data.mat'],'file')
        error('Download from http://audition.ens.fr/adc/NoiseTools/DATA/figures_deCheveigne_Arzounian_2017/');
    end
    load ([datadir,'fig_9_data']);
end

%% Add an artificial glitch on one channel:
% Subtract the mean of each channels:
x=nt_demean(x);
% 1000 mV amplitude glitch on channel 1 at 2 seconds:
x(1000+(1:100),1)=1000;

%% Preprocess signals:

% Detrend signals
[y,w]=nt_detrend(x,3,[],[],3); % get weight vector for robust fit of 3rd order polynomial
[y,w]=nt_detrend(x,10,w,[],3); % subtract 10th order polynomial

% Detect outliers automatically and replace them by signal interpolation:
[ww,yy]=nt_outliers(y,[],1);

%% Plot figure

% Prepare figure:
figure(1); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 150])

% Plot raw signal:
ch=79; % channel
plot(t, x(:,ch)+3);  hold on % raw signal
% Plot signal after standard re-referencing to the mean of all channels
plot(t,y(:,ch)-mean(y,2), 'k');
% Plot signal after robust re-referencing, i.e. using inpainted data:
plot(t,yy(:,ch)-mean(yy,2)-3,'r','linewidth',2);
% Set axes properties and legend:
ylim([-10 10])
set(gca,'fontsize',14, 'xtick', 1:4); ylabel('amplitude (mV)'); xlabel('time (s)');
legend('raw', 'rereferenced', 'robust r.', 'location', 'eastoutside'); legend boxoff

% Save figure
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_9')
