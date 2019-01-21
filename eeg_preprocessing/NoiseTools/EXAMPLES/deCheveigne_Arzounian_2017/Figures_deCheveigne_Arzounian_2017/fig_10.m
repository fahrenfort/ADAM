% Script to reproduce Fig. 10 of:
% de Cheveigné, A., Arzounian, D., Robust detrending, rereferencing,
% outlier detection, and inpainting for multichannel data,(submitted to
% NeuroImage).
%
% This script requires the NoiseTools toolbox
% (http://audition.ens.fr/adc/NoiseTools/). Tested with version
% 29-Nov-2017.

%% Load MEG data from file
clear;
datadir='../../../DATA/figures_deCheveigne_Arzounian_2017/';

if 0 % Create reduced MAT data file from raw source file   
    % Load ex1 ex2 ex3 stim:
    load alain3.mat 
    % Unfold MEG and stimulus signals:
    ex1=ex1(:);
    ex2=ex2(:);
    ex3=ex3(:);
    stim=stim(:);
    % Provide sampling rate:
    sr=2400; % (Hz)
    % Concatenate first 2 channels in cell array:
    ex={ex1,ex2};
    % Save relevant variables
    save ([datadir,'fig_10_data'], 'ex', 'stim', 'sr');
else % Load reduced MAT file
    if 2~=exist([datadir,'fig_10_data.mat'],'file')
        error('Download from http://audition.ens.fr/adc/NoiseTools/DATA/figures_deCheveigne_Arzounian_2017/');
    end
    load ([datadir, 'fig_10_data']);
end

%% Get the positions of stimulus onsets:
% Find high values:
stim=stim>0.5*10^5; 
% Pick first sample of each event:
stimIdx=find(~stim(1:end-1) & stim(2:end));
% Shift a bit to catch onset of step response:
ADVANCE=3; 
stimIdx=stimIdx-ADVANCE; 

%% Remove steps

% Set algorithm parameters:
GUARD=400; % minimum duration of a step
DEPTH=30; % determines how many steps can be removed
THRESH=0.1; % threshold variance reduction
MINSTEP=10; % minimum absolute step size 

% Take first channel:
x=ex{1};
% Convert from pT to fT:
x=x/1000;
% Remove steps:
xx=nt_destep(x,THRESH,GUARD,DEPTH,MINSTEP);

% Restrict signal to a few seconds:
FOCUS=201:4100; % target samples
t=(0:numel(FOCUS)-1)/sr; % Time vector
x=x(FOCUS);xx=xx(FOCUS); % reduced raw and de-stepped data vectors

% Subtract mean of signal:
x=nt_demean(x,1:100); % in raw signal
xx=nt_demean(xx,1:100); % in de-stepped signal

%% Plot results
% Prepare figure:
figure(1); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 300])
% Plot raw and de-stepped signals on the left:
subplot 121;
plot(t,x,'k','linewidth',0.5); hold on; % raw
plot(t,xx,'r','linewidth',1); hold on; % de-stepped
% Set legend and axes properties:
legend('raw','clean','location','northwest'); 
title('step removal'); legend boxoff 
set(gca,'fontsize',14,'ytick',[0 5000]); xlim([t(1) t(end)])
ylabel('field (fT)'); xlabel('time (s)');

%% Illustrate removal of ringing artifact:
% Take second channel:
x=ex{2};
% Convert from pT to fT:
x=x/1000;
% Remove steps:
xx=nt_destep(x,THRESH,GUARD,DEPTH,MINSTEP);
% Fit and remove ringing artifact:
xxx=nt_deboing(xx,stimIdx);
% Restrict signal to a few ms:
FOCUS=801:1200;
t=(0:numel(FOCUS)-1)/sr; % time vector
x=x(FOCUS);xx=xx(FOCUS);xxx=xxx(FOCUS); % reduced signals
% Subtract mean of signal:
x=nt_demean(x,1:100); % in raw signal
xx=nt_demean(xx,1:100); % in de-stepped signal
xxx=nt_demean(xxx,1:100); % in artifact-free signal
% Plot signals on the right of the figure:
subplot 122
plot(t,[xx],'k','linewidth',0.5); hold on; % de-stepped signal
plot(t,[xxx],'r','linewidth',1); hold on; % artifact-free
% Set axes properties and legend:
title('ringing removal'); 
legend('raw','clean','location','southeast'); 
legend boxoff
set(gca,'fontsize',14,'ytick',[-40 -20 0 20 40]); xlim([t(1) t(end)])
ylabel('field (fT)'); xlabel('time (s)')

%% Save 2 panels:
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_10');
