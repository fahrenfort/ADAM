% Script to reproduce Fig. 1 of:
% de Cheveigné, A., Arzounian, D., Robust detrending, rereferencing,
% outlier detection, and inpainting for multichannel data,(submitted to
% NeuroImage).
%
% This script requires the NoiseTools toolbox
% (http://audition.ens.fr/adc/NoiseTools/). Tested with version
% 29-Nov-2017.

clear
datadir='../../../DATA/figures_deCheveigne_Arzounian_2017/';

%% Load data from file
if 0 % Initial creation of MAT data file from raw source file
    
    % Ensure that biosig works
    if 2 ~=exist('sopen')
        error('download BioSig from http://biosig.sourceforge.net');
    end
    x = fileparts( which('sopen') );
    rmpath(x);
    addpath(x,'-begin');
    
    % Get data from raw source file and save
    fname='/data/meg/arzounian/ADC_DA_140521_p20/ADC_DA_140521_p20_01_calib';
    [p,x]=nt_read_data(fname);
    x=x'; % Reshape to time x channels
    DSR=8;
    x=nt_dsample(x,DSR);
    p.sr=p.sr/DSR;
    save ([datadir,'fig_1_data'], 'p', 'x');
    return
    
else % Load data from MAT file
    if 2~=exist([datadir,'fig_1_data.mat'])
        error('Download from http://audition.ens.fr/adc/NoiseTools/DATA/figures_deCheveigne_Arzounian_2017/');
    end
    load ([datadir,'fig_1_data.mat'])
end


%% Prepare data
x=nt_demean(x); % Remove mean of each channel
x=x/1000; % Convert from V to mV

%% Create figure
figure(2); clf; set(gcf,'color', [1 1 1], 'position',[100 100 395 213])
t=(0:size(x,1)-1)/p.sr; % time vector
plot(t,x);
set(gca,'fontsize',14);
xlabel('time (s)'); ylabel('amplitude (mV)'); 
xlim([t(1),t(end)]);
ylim([-4 4])

% Save figure
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_1')
