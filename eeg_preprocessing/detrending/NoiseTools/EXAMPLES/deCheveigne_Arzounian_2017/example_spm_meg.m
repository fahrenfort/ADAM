%{
This is an example script to illustrate routines of the NoiseTools toolbox.
- MEG data are loaded, downsampled, and the mean removed from each channel.
- Step glitches are removed, as well as the antialiasing filter step response. 
- Slow components and line power & harmonics are removed (3 different options).
- Eyeblink components are isolated and projected out.
- Data are cut into epochs based on triggers.
- DSS is applied to find repeatable components.
- A subset of components is selected and projected back to denoise the data.

Tested with NoiseTools version 03-Dec-2017.

Requires FieldTrip (or equivalent), and data that can be downloaded from SPM
website.
%}
clear

datadir='/data/meg/spm/MEG/';
savedir=[];%'./PP/'; % directory to save results

fname='SPM_CTF_MEG_example_faces1_3D.ds';
if ~exist([datadir,fname]);
    error('download data from http://www.fil.ion.ucl.ac.uk/spm/download/data/mmfaces/multimodal_meg.zip');
end

% read data
if ~2==exist('ft_read_header'); 
    disp('download FieldTrip from http://www.fieldtriptoolbox.org/download,'); 
    error('run ft_defaults to set path and add utilities/ directory');
end
%h=ft_read_header([datadir,fname]); % works on some systems, fails on others
h=ft_read_header([datadir,fname],'headerformat','ctf_old');
sr=h.Fs;
%x=ft_read_data([datadir,fname]); % works on some systems, fails on others
x=ft_read_data([datadir,fname],'headerformat','ctf_old','dataformat','ctf_old');
x=x';
trig=x(:,2);
x=x(:,4:306); % keep only MEG
ref=x(:,1:32);
x=x(:,33:end);

x=nt_demean(x);
ref=nt_demean(ref);

figure(1); clf; 
plot(x); title('step & ringing artifact removed'); drawnow

% remove step glitches (tricky because of linear trends)
thresh=0.8; % help nt_destep for explanation
guard=1000;
depth=6;
minstep=[];
for iChan=1:size(x,2);
    [y,stepList]=nt_destep(x(:,iChan),thresh,guard,depth,minstep);
    if ~isempty(stepList)
        stepList=[1,stepList,size(x,1)];
        for iStep=1:numel(stepList)-1
            x(stepList(iStep):stepList(iStep+1),iChan)=...
                nt_detrend(x(stepList(iStep):stepList(iStep+1),iChan),1);
        end
        tmp=nt_destep(x(:,iChan),thresh,guard,depth,minstep);
        % remove step response of antialiasing filter:
        ADVANCE=3; % shift step position forward to catch full IR
        x(:,iChan)=nt_deboing(tmp,stepList(iStep)-ADVANCE); clear tmp 
    end
end

x=x(1:end-5000,:); % clip glitch at very end
ref=ref(1:end-5000,:);

figure(2); clf;
plot(x); title('raw, demeaned'); drawnow

%{
There are at least three strategies to remove power line artifact and slow drifts:
(1) Smooth with a boxcar kernel of size sr/50Hz to remove 50Hz and harmonics. 
This is straightforward but the data are low-pass filtered (-3dB at 25 Hz).
The waveform is minimally affected (compared to standard low-pass or
notch). Apply polynomial detrending to remove slow drifts (as in EEG).
(2) Regress out the reference channels (this removes both 50 Hz and slow drifts).
Downside is that ref channels might pick up some brain activity (unlikely
if PCA is applied and only few PCs are selected).
(3) Use DSS to design a spatial filter to remove 50 Hz and slow components.
Use DSS to design a second spatial filter to remove very slow components.
%}

STRATEGY=1;
switch STRATEGY
    case 1
        figure(101); clf; 
        nt_spect_plot(x(1:10000,:),1024,[],[],sr); 
        x=nt_smooth(x,sr/50); % removes 50Hz and harmonics
        hold on; nt_spect_plot(x(1:10000,:),1024,[],[],sr); set(gca,'xscale','log')
        title('smoothing to suppress powerline'); legend('before', 'boxcar sr/50Hz'); drawnow
        ORDER=10;
        x=nt_detrend(x,ORDER); % remove slow artifact
    case 2
        z=nt_pca(ref); % PCA to reduce dimensionality
        NREMOVE=6;
        figure(101); plot(z(:,1:NREMOVE)); title ('PCs to remove');
        figure(102); clf; 
        nt_spect_plot(x(1:10000,:),1024,[],[],sr); 
        SHIFTS=1:10;
        x=nt_tsr(x,z(:,1:NREMOVE),SHIFTS);
        hold on; nt_spect_plot(x(1:10000,:),1024,[],[],sr); ; set(gca,'xscale','log')
        title('time-shift regression of ref channels'); legend('before', 'after'); drawnow
    case 3
        y=[ref,x];
        
        % remove slow components
        c0=nt_cov(y);
        c1=nt_cov(diff(y)); % emphasizes high frequencies, slow components will be last
        [todss,pwr0,pwr1]=nt_dss0(c0,c1);
        figure(101); clf; plot(pwr1./pwr0,'.-'); ylabel('score'); xlabel('components'); title('DSS slow components');
        z=nt_mmat(y,todss);
        z=fliplr(z); % put slow first
        NREMOVE=8; % based on examination of the cross-covariance matrix between z and ref        
        figure(102); plot(nt_normcol(z(:,1:NREMOVE))); title('slow components');
        y=nt_tsr(y,z(:,1:NREMOVE));
        x=nt_tsr(x,z(:,1:NREMOVE));
        
         
        % remove power line components
        c0=nt_cov(y);
        c1=nt_cov(nt_sparse_filter(y,[0; sr/50],[1; -1])); % bias suppresses DC & 50 Hz & harmonics
        [todss,pwr0,pwr1]=nt_dss0(c0,c1);
        figure(103); clf; plot(pwr1./pwr0,'.-'); ylabel('score'); xlabel('components'); title('DSS 50 Hz');
        figure(104); clf
        nt_spect_plot(x(1:10000,:),1024,[],[],sr); 
        z=nt_mmat(y,todss);
        z=fliplr(z); % put 50 Hz first       
        NREMOVE2=5;
        x=nt_tsr(x,z(:,1:NREMOVE2));
        hold on; nt_spect_plot(x(1:10000,:),1024,[],[],sr); set(gca,'xscale','log');
        title('DSS to remove 50 Hz components'); legend ('before', 'after');
end
        
figure(3); clf;
plot(x); title('after slow artifact & powerline removal'); drawnow


% remove eyeblink components
w=abs(x)> 4*median(abs(x(:)));
w=mean(w,2);
c0=nt_cov(x);
c1=nt_cov(bsxfun(@times,x,w));
[todss,pwr0,pwr1]=nt_dss0(c0,c1);
figure(4); clf; 
subplot 121; plot(pwr1./pwr0,'.-'); ylabel('score'); xlabel('components'); title('DSS eyeblink');
z=nt_mmat(x,todss);

% repeat to tighten estimate
w=abs(z(:,1))> 4*median(abs(z(:,1)));
c0=nt_cov(x);
c1=nt_cov(bsxfun(@times,x,w));
[todss,pwr0,pwr1]=nt_dss0(c0,c1);
subplot 122; plot(pwr1./pwr0,'.-'); ylabel('score'); xlabel('components'); title('DSS eyeblink 2');
z=nt_mmat(x,todss);

NREMOVE3=2;
x=nt_tsr(x,z(:,1:NREMOVE3));
figure(5); clf;
plot(x); title('after eyeblink removal'); drawnow

% cut into epochs
idxTrig=find(trig(1:end-1)<0.5 & trig(2:end)>0.5);
PRE=round(sr*0.3); % pad to anchor detrending
POST=round(sr*0.9); 
nsamples=PRE+POST;
nchans=size(x,2);
ntrials=numel(idxTrig);
xx=zeros(nsamples,nchans,ntrials); 
for iTrial=6:ntrials-5
    xx(:,:,iTrial)=x((idxTrig(iTrial)-PRE+1):(idxTrig(iTrial)+POST),:);
end
xx=nt_demean(xx);

% DSS to emphasize repeatablity
[todss,pwr0,pwr1]=nt_dss1(xx);
fromdss=pinv(todss);
z=nt_mmat(xx,todss);

figure(6); clf;
plot(pwr1./pwr0,'.-'); xlabel('component'); ylabel('score'); title ('repeatability DSS');

figure(7); clf;
t=linspace(-PRE/sr,POST/sr,nsamples);
nt_bsplot(z(:,1,:),[],[],t); title('best DSS component');
xlabel('s')

figure(8); clf
plot(t,nt_normcol(mean(z(:,1:7,:),3))); title('DSS components 1:7'); xlabel('s');

% denoise by selecting best components and projecting back to sensor space
NKEEP=7;
yy=nt_mmat(xx,todss(:,1:NKEEP)*fromdss(1:NKEEP,:));

figure(9); clf;
subplot 121; plot(t,mean(xx,3)); title('before DSS denoising'); xlabel('s')
subplot 122; plot(t,mean(yy,3)); title('after'); xlabel('s')


if ~isempty(savedir);
    if ~7==exist(savedir); mkdir(savedir); end
    [PATHSTR,NAME,EXT] = fileparts(fname);
    savename=[savedir,'/',NAME, '_PP'];
    save (savename, 'xx', 'sr', 'PRE', 'POST');
    eval(['!cp ', mfilename, '.m ',savedir]);  % save this script to document data
end


