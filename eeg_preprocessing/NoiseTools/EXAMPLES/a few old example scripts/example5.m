% This example uses DSS to isolate narrowband power near 16 Hz in MEG data.
% The presence of strong narrowband power near 10Hz would mask this
% component without DSS.
%
% The best component is dominated by 16 Hz, but the best sensor (ie the sensor 
% most strongly containing this component) is dominated mainly by 10 Hz.
%
% Uses nt_bias_filter(), nt_dss0().

clear;
disp(mfilename);
help(mfilename)

% load data
FNAME=[fileparts(which('nt_version')), '/DATA/example_data.mat'];
if ~exist(FNAME); 
    disp('file ''../DATA/example_data.mat'' not found, get it at http://cognition.ens.fr/Audition/adc/NoiseTools/DATA/');
    return
end
load(FNAME)  % loads 'meg', 'sr'
% excerpt of data from 
% Duncan, K.K., Hadjipapas, A., Li, S., Kourtzi, Z., Bagshaw, A., Barnes, G., 2009. Identifying
% spatially overlapping local cortical networks with MEG. Hum. Brain Mapp. 7,
% 1003?1016.
% With thanks to authors.

% first remove 50 Hz & harmonics (see example4)
disp('remove 50 Hz & harmonics...');
[c0,c1]=nt_bias_fft(meg,[50, 100, 150]/sr, 512);
[todss,pwr0,pwr1]=nt_dss0(c0,c1); 
p1=pwr1./pwr0; % score, proportional to power ratio of 50Hz & harmonics to full band
z=nt_mmat(meg,todss);
NREMOVE=20;
meg=nt_tsr(meg,z(:,1:NREMOVE,:)); % regress out to get clean data

% downsample
DSRATIO=3;
meg=nt_dsample(meg,DSRATIO);
sr=sr/DSRATIO;



% bias filter is second-order resonator
disp('DSS to isolate 16 Hz components...');
FPEAK=16; % Hz
Q=8; % determines width
[b,a]=nt_filter_peak(FPEAK/(sr/2),Q);

% covariance matrices of full band (c0) and filtered (c1)
[c0,c1]=nt_bias_filter(meg,b,a);

% DSS matrix
[todss,pwr0,pwr1]=nt_dss0(c0,c1); 
p1=pwr1./pwr0; % score, proportional to power ratio of 50Hz & harmonics to full band

% DSS components
z=nt_mmat(meg,todss);


% plot bias score
figure(1); clf; set(gcf,'color', [1 1 1]);
plot(p1, '.-'); xlabel('component'); ylabel('score'); title('DSS score');


% plot spectra of DSS components
figure(2); clf; set(gcf,'color', [1 1 1]);
nt_spect_plot2(nt_normcol(z(:,1:30,:)),512,[],[],sr);
title('spectra of first 30 DSS components'); ylabel('component')

% plot spectra of best DSS component and best sensor
figure(3); clf; set(gcf,'color', [1 1 1]);
nt_spect_plot(nt_normcol(z(:,1,:)),512,[],[],sr);
hold on;
[~,idx]=max(abs(nt_xcov(z(:,1,:),nt_normcol(meg))));
nt_spect_plot(nt_normcol(meg(:,idx,:)),512,[],[],sr);
nt_linecolors([], [3 1]);
legend('best component', 'best sensor'); legend boxoff; axis tight
title('bias to 16 Hz')


