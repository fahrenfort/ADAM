% Script to reproduce Fig. 3 of:
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
    fname='/data/meg/ADC_DA_140521_p20/ADC_DA_140521_p20_01_calib';
    [p,x]=nt_read_data(fname);
    % Reshape to time x channel:
    x=x'; 
    % Select electrode and time interval of interest:
    k=22; % electrode
    x=nt_demean(x(1800001:end-10000,k)); 
    % Save:
    save ([datadir,'fig_3_data'], 'p', 'x');
    
else % Load data from MAT file
    if 2~=exist([datadir, 'fig_3_data.mat'],'file')
        error('Download from http://audition.ens.fr/adc/NoiseTools/DATA/figures_deCheveigne_Arzounian_2017/');
    end
    load ([datadir, 'fig_3_data']);
end

%% Detrend the real signal
% Subtract mean from signal:
x=nt_demean(x); 
x=x/1000; % Convert from V to mV
% Fit polynomial and remove from signal:
ORDER=10; % Order of the polynomial
y=nt_detrend(x,ORDER);  % returns the detrended signal

% Plot results
figure(1); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 200])
t=(0:size(x,1)-1)/p.sr; % time vector
subplot 121;  % Plot original signal and trend on the left
% The trend is the difference between the original and detrended signals:
plot(t,x-y,'r','linewidth',2); hold on
plot(t,x, 'k'); 
set(gca,'fontsize',14,'ytick',[-.2 0 .2]); 
ylim([-0.25 0.25])
xlabel('time (s)'); xlim([t(1),t(end)])
ylabel('amplitude (mV)');
legend('fit','raw'); legend boxoff
subplot 122; % Plot detrended signal on the right
plot(t,y, 'k'); 
ylim([-0.25 0.25])
set(gca,'fontsize',14, 'ytick', []); 
xlabel('time (s)'); xlim([t(1),t(end)])
legend('detrended'); legend boxoff

% Save panels:
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_3_top')

%%  Add glitch and perform standard (i.e. non-robust) detrending
x(300000:400000)=0.8; % artificial glitch
% 'Non-robust' detrending is obtained by incorporating outlier samples in
% the fit:
thresh=100; % threshold for outliers, purposefully high here
y=nt_detrend(x,ORDER,[],[],thresh); % detrended signal

% Plot results
figure(2); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 200])
subplot 121;  % Plot original signal and trend on the left
plot(t,x-y,'r','linewidth',2); hold on
plot(t,x, 'k'); 
set(gca,'fontsize',14); 
ylim([-.5 1])
xlabel('time (s)'); xlim([t(1),t(end)])
ylabel('amplitude (mV)');
legend('fit','raw'); legend boxoff
subplot 122;  % Plot detrended signal on the right
plot(t,y, 'k'); hold on
plot(t,0*t,':'); % Add zero baseline
ylim([-.5 1])
set(gca,'fontsize',14, 'ytick', []); 
xlabel('time (s)'); xlim([t(1),t(end)])
legend('detrended'); legend boxoff

% Save panels
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_3_middle')

%% Perform robust detrending
% ... by discarding outlier samples from the fit:
w=ones(size(t')); % weighting function
w(300000:400000)=0; % Assign zero weight to glitch part
y=nt_detrend(x,ORDER,w); % detrended signal

% Plot results
figure(3); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 200])
subplot 121;  % Plot original signal and trend on the left
plot(t,x-y,'r','linewidth',2); hold on
plot(t,x, 'k'); 
% Display the weighting function underneath:
w(w==0)=nan;
plot(t,-0.4*w,'color',0.6*[1 1 1],'linewidth', 4);
set(gca,'fontsize',14); 
ylim([-.5 1])
xlabel('time (s)'); xlim([t(1),t(end)])
ylabel('amplitude (mV)');
legend('fit (robust)','raw', 'mask'); legend boxoff
subplot 122; % Plot detrended signal on the right
plot(t,y, 'k'); hold on
plot(t,0*t,':'); % Add zero baseline
ylim([-.5 1])
set(gca,'fontsize',14, 'ytick', []); 
xlabel('time (s)'); xlim([t(1),t(end)])
legend('detrended (robust)'); legend boxoff

% Save panels
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_3_bottom')


    