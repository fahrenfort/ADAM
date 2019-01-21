% Script to reproduce Fig. 6 and 7 of:
% de Cheveigné, A., Arzounian, D., Robust detrending, rereferencing,
% outlier detection, and inpainting for multichannel data,(submitted to
% NeuroImage).
%
% This script requires the NoiseTools toolbox
% (http://audition.ens.fr/adc/NoiseTools/). Tested with version
% 29-Nov-2017.

clear

% [Need to set seed  for random number generation ?]

%% Create artificial, corrupted signal mixture:
nsamples=1000; % number of samples
nchans=50; % number of channels
N=10; % number of sources
t=(0:nsamples-1)/100; % time vector
% Create source signals:
x=zeros(nsamples,N); 
for k=1:N
    x(:,k)=sin(2*pi*k*(1:nsamples)'/nsamples);
end
% Mix signals using random, linear projection matrix:
x=x*randn(N,nchans);
% Make a few channels stronger:
x(:,1:6)=x(:,1:6)*3;
% Store a copy of the clean channel signals:
x0=x;
% Add glitches:
w=ones(size(x)); % Initialize the weight matrix (all 1 means all samples allright so far);
glitchwidth=20; % glitch duration in samples
for k=1:nchans % browse channels one by one
    % Select a random position in time:
    idx=100+ceil(rand*(nsamples-300)); 
    % Add an offset and random noise to the signal at this location:
    x(idx+(1:glitchwidth),k)=x(idx+(1:glitchwidth),k)+40+10*randn(glitchwidth,1);  
    % Update weight matrix to signal the glitch:
    w(idx+(1:glitchwidth),k)=0; 
end

%% Use the inpainting algorithm to replace glitches with an interpolated signal:
y=nt_inpaint(x,w);

%% Plot figure

% Prepare figure
figure(1); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 400])

% Plot raw signals at the top:
subplot 311;
plot(t,x, 'linewidth',0.3); hold on; % Corrupted signal
plot(t,x0, 'linewidth',1); title ('raw') % Clean signal on top
% Set axes properties:
set(gca,'fontsize',12,'ytick',[],'xticklabel', []); ylabel('A.U.') 
ylim([min(x(:)),max(x(:))]*1.3)

% Plot weight matrix in the middle:
subplot 312;
nt_imagescc(2*w'-1);
% Set axes and colormap properties:
set(gca,'fontsize',12,'ytick',[1 50 100],'xtick', [],'xticklabel', [0 1 2 3]); %xlabel('time (s)');
ylabel('channel')
c=colormap; colormap([repmat(c(1,:),32,1);repmat(c(end,:),32,1)]); title('weight (given)')
% Add colorbar on the right side and set properties:
h=colorbar('location','east'); 
set(h,'limits',[-1 1], 'ytick',[-1 1],'yticklabel',[0 1],'fontsize',12)

% Plot inpainted signals at the bottom:
subplot 313;
plot(t,x0*nan, 'linewidth',1); hold on % to get same colors as first plot
plot(t,y, 'linewidth',1);
% Set axes properties:
set(gca,'fontsize',12,'ytick',[],'xtick', 0:15); xlabel('time (s)'); ylabel('A.U.'); title ('interpolated')
ylim([min(x(:)),max(x(:))]*1.3)

% Save 3 panels:
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_6')

%% Perform inpainting of corrupted signals, but detecting glitches automatically:
[ww,yy]=nt_outliers(x,[],1,10);

%% Plot new results in new figure

% Prepare figure:
figure(2); clf; set(gcf,'color', [1 1 1], 'position',[100 100 600 270])

% Plot the estimated weight matrix at the top:
subplot 211;
nt_imagescc(2*ww'-1);
% Set colormap and axes properties as in previous figure
set(gca,'fontsize',12,'ytick',[1 50 100],'xtick', [],'xticklabel', [0 1 2 3]); %xlabel('time (s)');
ylabel('channel')
c=colormap; colormap([repmat(c(1,:),32,1);repmat(c(end,:),32,1)]); title('weight (estimated)')
% Add colorbar
h=colorbar('location','east');
set(h,'limits',[-1 1], 'ytick',[-1 1],'yticklabel',[0 1],'fontsize',12)

% Plot inpainted signals at the bottom:
subplot 212;
plot(t,yy, 'linewidth',0.3); hold on; % to get same colors as first plot
plot(t,yy, 'linewidth',1); title ('interpolated (using estimated weight)')
% Set axes properties:
set(gca,'fontsize',12,'ytick', [], 'xtick', 0:15); xlabel('time (s)'); ylabel('A.U.')
ylim([min(x(:)),max(x(:))]*1.3)

% Save two panels:
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', 'Fig_7')


