%{
This is an example script to illustrate routines of the NoiseTools toolbox.
- EEG data are loaded, downsampled, and the mean removed from each channel.
- Bad channels are detected and interpolated based on neighbors.
- The slow varying trend is removed by robust fitting with a polynomial.
- Optionally: corrupt samples are detected and replaced based on intact data.
- Eyeblink components are isolated and projected out.
- Data are rereferenced by subtracting the robust mean.
- Data are cut into epochs based on triggers.
- The residual trend is removed by robust fitting each trial & channel with a linear function.
- Outlier trials are detected and discarded.
- DSS is applied to find repeatable components.
- A subset of components is selected and projected back to denoise the data.

Tested with NoiseTools version 03-Dec-2017.

Requires FieldTrip (or equivalent), and data that can be downloaded from SPM
website.
%}
clear

datadir='/data/meg/spm/EEG/';
savedir=[];%'./PP/'; % directory to save results

fname='faces_run2.bdf';
if ~exist([datadir,fname]);
    error('download data from http://www.fil.ion.ucl.ac.uk/spm/download/data/mmfaces/multimodal_eeg.zip');
end

% read data
h=ft_read_header([datadir,fname]); 
x=ft_read_data([datadir,fname]);
x=x';
trig=x(:,144);
x=x(:,1:128);
sr=h.Fs;

% [p,x]=nt_read_data([datadir,fname]);
% trig=x(:,144);
% x=x(:,1:128);
% sr=p.sr;

% downsample to save computation
DSR=8; % --> 256 Hz
x=nt_dsample(x,DSR);
x=nt_demean(x);
sr=sr/DSR;

figure(1); clf
plot(x); title('raw, demeaned');

% read electrode position info file
fid = fopen([datadir,'electrode_locations_and_headshape.sfp'], 'r');
coordinates=textscan(fid, '%s %f %f %f ');
coordinates=cell2mat(coordinates(2:end));
coordinates=coordinates(4:131,:);
fclose(fid);
proximity_matrix=nt_proximity(coordinates); % between electrodes

% find bad channels, interpolate
proportion=0.5; % criterion proportion of bad samples
thresh1=3; % threshold in units of median absolute value over all data
thresh2=[]; % absolute threshold 
thresh3=[]; % absolute threshold applied to sns-processed data
iBad=nt_find_bad_channels(x(:,1:128),proportion,thresh1,thresh2,thresh3);
[toGood,fromGood]=nt_interpolate_bad_channels(x(:,1:128),iBad,proximity_matrix);
x(:,1:128)=x(:,1:128)*(toGood*fromGood);
% iBad=nt_find_bad_channels(x(:,129:end),proportion,thresh1,thresh2,thresh3);
% x(:,128+iBad)=0; % positions unknown, just zap


figure(2); clf
plot(x); title('bad channels interpolated');

% detrend
ORDER=10; % of polynomial
[x,w]=nt_detrend(x,ORDER);
figure(3); clf
plot(x); title('detrended');
figure(4); clf
nt_imagescc(w'); title('weights from nt_detrend','interpreter','none'); ylabel('channels')


% find outliers, inpaint
if 0 % not needed for these data, time consuming
    
    thresh=2; % threshold for declaring an outlier
    niter=3; % number of iterations
    [w,y]=nt_outliers(x,w,thresh,niter);
    
    
    figure(101); clf;
    nt_imagescc(w'); title('weights from nt_outliers'); ylabel('channels')
    
    if 1
        x=y; % all samples replaced
    else
        y=nt_inpaint(x,w); % interpolate over outliers
        x=y; % only outlier samples replaced
    end
    
    figure(102); clf
    plot([x,y]); legend('raw','outliers interpolated');
end

if 1
    % remove eyeblinks
    eye_channels=[93 94 81 82 71 72];
    [B,A]=butter(2,1/(sr/2), 'high');
    tmp=nt_pca(filter(B,A,x(:,eye_channels)));
    mask=abs(tmp(:,1))>3*median(abs(tmp(:,1)));
    C0=nt_cov(x);
    C1=nt_cov(bsxfun(@times, x,mask));
    [todss,pwr0,pwr1]=nt_dss0(C0,C1);
    figure(5); plot(pwr1./pwr0, '.-'); ylabel('score'); xlabel('component'); title ('eyeblink DSS');
    eye_components=x*todss;
    NREMOVE=2; % arbitrary, may need adjusting
    x=nt_tsr(x,eye_components(:,1:NREMOVE));
    figure(6); clf;
    subplot 211; plot(eye_components(:,1:NREMOVE)); title('eye components');
    subplot 212; plot(x); title('eyeblinks removed');
end

% rereference
x=nt_rereference(x,w);


% cut into epochs
trig=nt_demean(trig);
trig(1:5000)=0; % zap initial pulse
trig=nt_demean(trig);
idxTrig=find(trig(1:end-1)<0.5 & trig(2:end)>0.5);
idxTrig=round(idxTrig/DSR);
PRE=round(sr*0.5); % pad to anchor detrending
POST=round(sr*1); 
nsamples=PRE+POST;
nchans=size(x,2);
ntrials=numel(idxTrig);

xx=zeros(nsamples,nchans,ntrials); 
for iTrial=6:ntrials-5
    xx(:,:,iTrial)=x((idxTrig(iTrial)-PRE+1):(idxTrig(iTrial)+POST),:);
end

% remove linear trend from each epoch (not really needed here)
xx=nt_detrend(xx,1);
xx=xx(PRE+1:PRE+1 + round(0.6*sr), :,:);  % chop off padding
xx=nt_demean(xx);

% remove outlier trials
THRESH=5;
iKeep=nt_find_outlier_trials(xx,THRESH);
xx=xx(:,:,iKeep);
xx=nt_demean(xx);

% DSS to emphasize repeatablity
[todss,pwr0,pwr1]=nt_dss1(xx);
fromdss=pinv(todss);
z=nt_mmat(xx,todss);

figure(9); clf;
subplot 121;
plot(pwr1./pwr0,'.-'); xlabel('component'); ylabel('score'); title ('repeatability DSS');
subplot 122;
nt_bsplot(z(:,1,:)); title('best DSS component');

% denoise by selecting best components and projecting back to sensor space
NKEEP=7;
y=nt_mmat(xx,todss(:,1:NKEEP)*fromdss(1:NKEEP,:));

figure(10); clf;
subplot 121; plot(mean(xx,3)); title('before DSS denoising');
subplot 122; plot(mean(y,3)); title('after');


if ~isempty(savedir)
    if ~7==exist(savedir); mkdir(savedir); end
    [PATHSTR,NAME,EXT] = fileparts(fname);
    savename=[savedir,'/',NAME, '_PP'];
    save (savename, 'xx', 'sr', 'PRE', 'POST');
    eval(['!cp ', mfilename, '.m ', savedir]); % save this script to document data
end

