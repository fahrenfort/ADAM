% Find linear combination of LFPs that maximizes selectivity to stimulus
% frequency.  Data are from a 16-channel electrode array in guinea pig
% auditory cortex. Stimuli are tone pips with frequencies from 0.56 to 36
% kHz in 6% steps (97 frequencies), each presented 8 times in pseudorandom
% order.
%
% Uses nt_dss0().

clear;
disp(mfilename);
help(mfilename)

FNAME=[fileparts(which('nt_version')), '/DATA/Fichier_1815_lfp.mat'];
if ~exist(FNAME); 
    disp(['file ',FNAME,' not found, get it at http://cognition.ens.fr/Audition/adc/NoiseTools/DATA/']);
    return
end

load (FNAME);   % loads LFPSampleRate, frequencies, xx
x=xx; clear xx
sr=LFPSampleRate;
[nsample,nchan,ntrial,nfreq]=size(x);
t=(0:nsample-1)'/sr;

% average over repeats, trim to first 100 ms
x=squeeze(mean(x,3)); % --> time * channels * frequency
x=x(find(t<=0.1),:,:,:);
nsample=size(x,1);
x=nt_demean(x);


% Component analysis:
% For each frequency, find the matrix that defines the linear combination of channels that
% maximizes response at that frequency relative to all others
to_otc=zeros(nchan,nfreq);
c0=nt_cov(x); 
for iBias=1:nfreq
    c1=nt_cov(x(:,:,iBias));
    [todss]=nt_dss0(c0,c1);
    todss=todss(:,1); % keep best
    to_otc(:,iBias)=todss;
end

% apply that matrix to get components
otc=zeros(nsample,nfreq,nfreq);    
for iBias=1:nfreq
    z=nt_mmat(x,to_otc(:,iBias));
    otc(:,iBias,:)=z; % optimally tuned component (time*bias*freq)
end

% for convenience, flip signs so all components are similar
a=nt_unfold(otc);
aa=nt_pca(a);
for iBias=1:nfreq
    if aa(:,1)'*a(:,iBias)<0
        otc(:,iBias,:)=-otc(:,iBias,:);
        to_otc(:,iBias)=-to_otc(:,iBias);
    end
end



% Display results.

% Tuning curves for electrodes and components

% RMS over time --> stim bias freq * stim freq
zz=squeeze(sqrt(mean(otc.^2,1)));
zz=nt_normcol(zz);

% same for electrode signals --> stim freq * electrodes
xx=squeeze(sqrt(mean(x.^2)))';
xx=nt_normcol(xx);

figure(1); clf
subplot 211
nt_imagescc(xx'); title('electrode tuning curves');
xlabel('stim freq (.56 - 36 kHz)');
ylabel('electrode');

subplot 212
nt_imagescc(zz); title('optimally tuned components');
xlabel('stim freq (.56 - 36 kHz)');
ylabel('bias freq (.56 - 36 kHz)');


figure(2); clf
plot([xx(:,11), zz(:,70)]);
xlabel('stim freq (.56 - 36 kHz)');
nt_linecolors([],[1 3]);
legend('electrode #11','component #70', 'location', 'northwest'); legend boxoff
title('example tuning curves');

figure(3); clf
subplot 121; 
nt_imagescc(squeeze(x(:,11,:)))
xlabel('stim freq (.56 - 36 kHz)');
ylabel('time (0-100ms)');
title('electrode 11');
subplot 122; 
nt_imagescc(squeeze(otc(:,72,:)))
xlabel('stim freq (.56 - 36 kHz)');
ylabel('time (0-100ms)');
title('component 72');

figure(4); clf
nt_imagescc(to_otc)
xlabel('bias freq (.56 - 36 kHz)');
ylabel('electrode'); 
title('weights for each component');

