% demo code for MCCA

disp('EEG data: basic MCCA');

%{
EEG data in response to repeated presentation of a 1 kHz tone for 25 subjects.

Data are organized as a 4D matrix with dimensions time, channel, trial, subject.

Data are preprocessed prior to epoching. Bad channels are detected automatically and removed, 
and data are downsampled to 128 Hz and filtered or detrended (several versions).

There are at least 124 repeats for each subject. We can use this to know
what is the 'true' cortical response to stimulation, and use it a 'ground
truth' to evaluate MCCA.

%}


clear
fname='../../DATA/MCCA_EXAMPLE_DATA/localiser_hpf_1Hz.mat';
if exist(fname, 'file')
    load (fname); % loads x sr preStim postStim 
else
    error('Download example data from http://audition.ens.fr/adc/NoiseTools/DATA/MCCA_EXAMPLE_DATA/');
end

[nsamples,nchans,ntrials,nsets]=size(x);
t=linspace(-preStim/sr, postStim/sr, nsamples);

% remove mean over pre-stimulus segment ('baseline correction')
for iSet=1:nsets
    x(:,:,:,iSet)=nt_demean(x(:,:,:,iSet),find(t<0));
end

% reduce dimensionality
NKEEP=20; 
y=zeros(nsamples,NKEEP,ntrials,nsets);
for iSet=1:nsets
    a=nt_pca(x(:,:,:,iSet)); 
    a=a(:,1:NKEEP,:);
    y(:,:,:,iSet)=a;
end

% unfold each data matrix (concatenate trials)
yy=permute(y,[1 3 2 4]);
yy=reshape(yy,[nsamples*ntrials,NKEEP,nsets]);

% concatenate over channels
yy=yy(:,:);

% MCCA
C=yy'*yy;
[A,score,AA]=nt_mcca(C,NKEEP);
z=yy*A;

% fold SCs (--> samples * SCs * trials)
z=nt_fold(z,nsamples);

figure(1);
subplot 121;
plot(score, '.-'); title('variance profile'); xlabel('SC'); ylabel('variance');
subplot 222;
nt_bsplot(z(:,1,:),[],[],t);
title('SC 1, trial avg'); set(gca,'ytick',[]); xlabel('time re stim onset (s)'); 
legend('+/- 2 SE', 'mean', 'location', 'northwest'); legend boxoff; 
subplot 224
nt_imagescc(nt_normcol(mean(z(:,1:30,:),3))');
title('first 30 SCs, trial avg');
ylabel('SC'); xlabel('time (samples)');
% Fig1: SCs found by MCCA are dominated by stimulus locked activity

%{
Trial-averaged SCs share quite similar time-courses, despite the fact that
non-averaged SCs are mutually orthogonal.  To summarize theseS time courses 
and remove redundancy we apply PCA to selected time-averaged SCs.  The PCA
matrix is then applied to the unaveraged SCs.
%}

% PCA selected trial-averaged SCs
NKEEP2=100;
topcs=nt_pca0(mean(z(:,1:NKEEP2,:),3)); 

% apply PCA matrix to unaveraged SCs
zz=nt_mmat(z(:,1:NKEEP2,:),topcs);

% tweak component signs (aesthetics)
sgn=sign(mean(mean(zz,3)));
sgn=reshape(sgn,[1 numel(sgn) 1]);
zz=bsxfun(@times,zz,sgn(1:NKEEP2));

figure(2);
subplot 121;
nt_imagescc(mean(zz(:,1:10,:),3)');
title('first 10 PCs of trial-avg SCs'); xlabel('time (samples)'); ylabel('PC');
for k=1:4;
    subplot (4,2,k*2);
    nt_bsplot(zz(:,k,:),[],[],t); 
    set(gca,'ytick',[]); title(['PC ',num2str(k)])
end
xlabel('time (s)');
%Fig2: The stimulus-locked activity of all subjects seems to be spanned by at least
% 4 components. These reflect both the potential multiplicity of
%stimulus-locked sources within each subject, and between-subject
%differences (e.g. in latency).



