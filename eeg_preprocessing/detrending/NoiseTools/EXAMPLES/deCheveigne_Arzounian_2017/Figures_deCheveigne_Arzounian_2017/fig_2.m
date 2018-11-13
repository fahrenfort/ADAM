% Script to reproduce Fig. 2 of:
% de Cheveigné, A., Arzounian, D., Robust detrending, rereferencing,
% outlier detection, and inpainting for multichannel data,(submitted to
% NeuroImage).
%
% This script requires the NoiseTools toolbox
% (http://audition.ens.fr/adc/NoiseTools/). Tested with version
% 29-Nov-2017.

clear

%% Create artificial pulse
x=zeros(400,1); % Zero baseline
x=[x;sin(pi*linspace(0,1,50)')]; % Pulse = half sinusoid
x=[x;zeros(400,1)]; % Return to baseline
sr=100; % Sampling rate (Hz)
t=(0:(size(x,1)-1))/sr; % Time vector
t=t-400/sr; % Shift 0 at pulse onset

%% Filter signal with different high-pass filters
orders={2,2,6}; % Filter orders
cutoffs={1,0.1, 0.1}; % Cutoff frequencies (Hz)
yy=[];yyy=[];
for k=1:numel(orders)
    % Set filter cutoff from frequency
    wn=cutoffs{k}/(sr/2); % cutoff
    % Design Butterworth high-pass filter with desired order
    [B,A]=butter(orders{k}, wn, 'high');
    % Filter signal causally (forward only)
    yy(:,k)=filter(B,A,x);
    % Filter signal non-causally (forward and backward)
    yyy(:,k)=filtfilt(B,A,x);
end

%% Plot results
figure(1); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 200])
axes('position', [12,20,30,75]/100)
% Plot true signal
plot(t,[x], 'k', 'linewidth', 2); hold on
% Plot outputs of causal filters underneath
plot(t,[yy(:,1)]-1, 'r', 'linewidth', 1);
plot(t,[yy(:,2)]-1, 'g', 'linewidth', 1);
plot(t,[yy(:,3)]-1, 'b', 'linewidth', 1);
% Plot outputs of noncausal filters underneath
plot(t,[yyy(:,1)]-2.7, 'r', 'linewidth', .5);
plot(t,[yyy(:,2)]-2.7, 'g', 'linewidth', .5);
plot(t,[yyy(:,3)]-2.7, 'b', 'linewidth', .5);
% Add zero reference lines
line([-10 10],[0 0], 'linestyle', ':')
line([-10 10],[-1 -1], 'linestyle', ':')
line([-10 10],[-2.7 -2.7], 'linestyle', ':')
line([0 0], [-5 5], 'linestyle', ':') 
% Label axes
set(gca,'fontsize',14, 'ytick', [], 'xtick',[-5:5]); 
ylim([-3.2, 1.3]); 
xlim([-1.2 2.2])
xlabel('time (s)'); ylabel('a.u.')

%% Create artificial step
x=zeros(20*sr,1); % Zero baseline
x=[x;0.5-0.5*cos(pi*linspace(0,1,50)')]; % Steep ramp
x=[x;ones(30*sr,1)]; % Plateau
t=linspace(-20, 30, size(x,1)); % New time vector

%% Filter signal with same high-pass filters as previously
yy=[];yyy=[];
for k=1:numel(orders)
    % Set filter cutoff from frequency
    wn=cutoffs{k}/(sr/2); % cutoff
    % Design Butterworth high-pass filter with desired order
    [B,A]=butter(orders{k}, wn, 'high');
    % Filter signal causally (forward only)
    yy(:,k)=filter(B,A,x);
    % Filter signal non-causally (forward and backward)
    yyy(:,k)=filtfilt(B,A,x);
end

%% Plot results in new sub-plot
axes('position', [47,20,45,75]/100)
% Plot true signal
plot(t,[x], 'k', 'linewidth', 2); hold on
% Plot outputs of causal filters underneath
plot(t,[yy(:,1)]-1, 'r', 'linewidth', 1);
plot(t,[yy(:,2)]-1, 'g', 'linewidth', 1);
plot(t,[yy(:,3)]-1, 'b', 'linewidth', 1);
% Plot outputs of noncausal filters underneath
plot(t,[yyy(:,1)]-2.5, 'r', 'linewidth', .5);
plot(t,[yyy(:,2)]-2.5, 'g', 'linewidth', .5);
plot(t,[yyy(:,3)]-2.5, 'b', 'linewidth', .5);
% Add zero reference lines
line([-50 50],[0 0], 'linestyle', ':')
line([-50 50],[-1 -1], 'linestyle', ':')
line([-50 50],[-2.5 -2.5], 'linestyle', ':')
line([0 0], [-5 5], 'linestyle', ':') 
% Label axes and add legend
set(gca,'fontsize',14, 'ytick', [], 'xtick',[-30:10:30]); 
ylim([-3.2, 1.3]); 
xlim([-12 22])
xlabel('time (s)'); ylabel('a.u.')
h=legend ('true signal','1 Hz, N=2','.1 Hz, N=2', '.1 Hz, N=6', 'location', 'eastoutside'); legend boxoff;

%% Save figure pannels
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_2_top')

%% Create artificial ramp
x=linspace(1,-1,50*sr)';
t=(0:(size(x,1)-1))/sr; % Time vector

%% Filter signal with same high-pass filters as previously
yy=[];yyy=[];
for k=1:numel(orders)
    % Set filter cutoff from frequency
    wn=cutoffs{k}/(sr/2); % cutoff
    % Design Butterworth high-pass filter with desired order
    [B,A]=butter(orders{k}, wn, 'high');
    % Filter signal causally (forward only)
    yy(:,k)=filter(B,A,x);
    % Filter signal non-causally (forward and backward)
    yyy(:,k)=filtfilt(B,A,x);
end

%% Plot results in new figure
figure(2); clf; set(gcf,'color', [1 1 1], 'position',[100 100 470 170])
% Plot true signal
plot(t,[x]+1, 'k', 'linewidth', 2); hold on
% Add zero reference line
line([-100 100],[1 1],'linestyle', ':');
% Plot outputs of causal filters
plot(t,[yy(:,1)]-1.5, 'r', 'linewidth', 1);
plot(t,[yy(:,2)]-1.5, 'g', 'linewidth', 1);
plot(t,[yy(:,3)]-1.5, 'b', 'linewidth', 1);
% Label axes
set(gca,'fontsize',14, 'ytick', []); 
ylim([-2.5 2.2]); 
xlim([-1 t(end)+1]); 
xlabel('time (s)'); ylabel('a.u.');
% Save figure panel
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_2_bottom')
