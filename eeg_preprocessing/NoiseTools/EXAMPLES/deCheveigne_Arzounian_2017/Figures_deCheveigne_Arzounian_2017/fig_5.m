% Script to reproduce Fig. 5 of:
% de Cheveigné, A., Arzounian, D., Robust detrending, rereferencing,
% outlier detection, and inpainting for multichannel data,(submitted to
% NeuroImage).
%
% This script requires the NoiseTools toolbox
% (http://audition.ens.fr/adc/NoiseTools/). Tested with version
% 29-Nov-2017.

%% Load EEG data from file for top panels
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
    fname='/DATA/MEG/LeGal/ADC_MLGDK_150619_pd01/ADC_MLGDK_150619_pd01_02.bdf';
    h=sopen(fname);
    x=sread(h);
    % Down-sample signal
    DSR=20; % Down-sampling ratio
    sr=h.SampleRate/DSR; % New sampling rate
    x=nt_dsample(x,DSR); % Down-sampled signal
    % Subtract the mean of each channel
    x=nt_demean(x);
    % Save relevant data
    save ([datadir, 'fig_5_data'], 'x', 'sr', 'fname');
else 
    load ([datadir, 'fig_5_data']);
end

%% Prepare signal and perform robust detrending:
% Keep a smaller time interval to zoom a bit:
x=nt_demean(x(1:21000,:));
% Convert from V to mV:
x=x/1000;
% Perform robust detrending on channel of interest:
ch=1; % Channel number
[y,w]=nt_detrend(x(:,ch),30,[],[],3.5); % Detrended signal
z=x(:,ch)-y; % Removed polynomial trend
t=(0:size(x,1)-1)'/sr; % Time vector

%% Plot raw and detrended signals:
% Prepare figure
figure(1); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 200])
 
% Plot raw signal and trend on the left:
subplot 121;
plot(t,z,'r','linewidth', 2); % trend
hold on; plot(t,[x(:,ch)], 'k'); % raw signal
% Set axes properties:
mx=max(x(:,ch)); mn=min(x(:,ch));
ylim([-.5 .5]); xlim([t(1)-3 t(end)+3]);
set(gca,'fontsize',14,'ytick',[-.5 0 .5]); xlabel('time (s)'); ylabel('amplitude (mV)');
title('raw')

% Plot detrended signal on the right:
subplot 122;
plot(t,nt_demean(y), 'k')
% Set axes proporties:
ylim([-.5 .5]); xlim([t(1)-3 t(end)+3]);
set(gca,'fontsize',14,'ytick',[]); xlabel('time (s)')
title('robust detrend')

% Save two panels:
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_5_top')

%% Create artificial, corrupted sinusoidal signal for bottom pannels:
clear
N=10000; % Number of samples
sr=10000; % Sampling rate (Hz)
t=(0:N-1)'/sr; % Time vector
% The noise is a 50-Hz sinusoid:
noise=sin(2*pi*50*t);
% The clean signal is a 1.25-Hz sinusoid
signal=-sin(2.5*pi*t);
% The signal is a weighted sum of clean signal and noise:
x=4*signal+noise;
% Add short glitch:
x(1001:1100)=200;

%% Plot raw signal and detrended versions:

% Prepare figure:
figure(2); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 213])

% Plot raw signal on the left:
axes('position', [12,20,25,70]/100)
plot(t,x, 'k');
hold on; plot(t,0*t,':k') % Add zero baseline
% Set axes properties:
set(gca,'fontsize',14,'ytick', [-10 0 10], 'xtick', [0 1]); title('raw'); xlabel('time (s)'); ylabel('A.U.')
ylim([-10 10]);

% Plot result of standard (i.e. non-robust) detrending in the middle:
axes('position', [42,20,25,70]/100)
% Detrending consists here in removing a 50-Hz sinusoidal trend. We
% therefore use the noise signal as the function basis. Non-robust
% detrending is obtained using a high threshold for outliers, e.g. 1000
% here:
plot(t,nt_detrend(x,[],[],noise,1000), 'k');
% Set axes properties:
hold on; plot(t,0*t,':k')
set(gca,'fontsize',14,'ytick',[], 'xtick', [0 1]); title('detrend'); xlabel('time (s)'); 
ylim([-10 10]);

% Plot result of robust detrending on the right:
axes('position', [72,20,25,70]/100)
% Same detrending as previouly, but using a lower threshold, e.g. 3 here:
plot(t,nt_detrend(x,[],[],noise,3), 'k');
% Set axes properties:
hold on; plot(t,0*t,':k')
set(gca,'fontsize',14,'ytick',[], 'xtick', [0 1]); title('robust detrend'); xlabel('time (s)'); 
ylim([-10 10]);

% Save 3 panels:
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_5_bottom')
