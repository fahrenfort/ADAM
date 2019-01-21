function [x,w,ww]=nt_star2(x,thresh,closest,w)
% [y,w,ww]=nt_star2(x,thresh,closest,w) - sensor noise suppression
%
%  y: denoised data 
%  w: 1 for parts used to calculate covariance (time*1)
%  ww: 0 for parts that needed fixing, 1 elsewhere (time*chans)
%
%  x: data to denoise (time*chans or time*chans*trials)
%  thresh: threshold for excentricity measure (default:1);
%  closest: indices of channels that are closest to each channel (default: all)
%  w: initial covariance weight matrix (time*1)
% 
% See also: nt_sns, nt_proximity

nt_greetings;

PCA_THRESH=10^-15;  % threshold for discarding weak PCs
NSMOOTH=10;         % samples, smoothing to apply to excentricity
MINPROP=0.4;        % minimum proportion of artifact free at first iteration
NITER=4;            % iterations to refine c0
VERBOSE=1;          % set to 0 to shut up
THRESH_RATIO=2;     % ratio of shared-weight and per-channel-weight thresholds

if nargin<1; error; end
if nargin<2 || isempty(thresh); thresh=1; end
if nargin<3; closest=[]; end
if ~isempty(closest)&&size(closest,1)~=size(x,2);
    error('closest array should have as many rows as channels of x'); 
end
if nargin<4; w=[]; end

if nargout==0 % plot, don't return result
    [y,w,ww]=nt_star2(x,thresh,closest,w);
    disp([mean(w(:)), mean(ww(:))])
    figure(1); clf;
    subplot 411; plot(x); mx=max(x(:)); mn=min(x(:)); ylim([mn mx]); title('raw');
    subplot 412; plot(y);ylim([mn,mx]); title('clean')
    subplot 413; plot(x-y); title('diff');
    subplot 414; plot(w, '.-'); ylim([0 1.1]); title('weight');
    clear x w ww
    return
end

[nsample,nchan,~]=size(x);
x=nt_unfold(x);

p0=nt_wpwr(x);
mn=mean(x); % save means
x=nt_demean(x);
nn=sqrt(mean(x.^2)); % save norm
x=nt_normcol(x);
p00=nt_wpwr(x);

% NaN, zero, and zero-weight channels are set to rand, which effectively excludes them
iNan=find(all(isnan(x)));
iZero=find(all(x==0));
x(:,iNan)=randn(size(x,1),numel(iNan));
x(:,iZero)=randn(size(x,1),numel(iZero));

x=nt_demean(x);

if isempty(w); w=ones(size(x,1),1); end
if size(w,2)>1; w=min(w,[],2); end
c0=nt_cov(x,[],w);          % initial covariance estimate
ww=ones(size(x));           % per-channel weights

%{
Find time intervals where at least one channel is excentric --> w==0.
%}

iIter=NITER;
while iIter>0
      
    for iChan=1:nchan

        % other channels
        if ~isempty(closest); 
            oChan=closest(iChan,:);
        else
            oChan=setdiff(1:nchan,iChan);
        end
        oChan(oChan>nchan)=[];
        
        % PCA other channels to remove weak dimensions
        [topcs,eigenvalues]=nt_pcarot(c0(oChan,oChan)); % PCA
        idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
        topcs=topcs(:,idx);
        
        % project this channel on other channels
        b=c0(iChan,oChan)*topcs/(topcs'*c0(oChan,oChan)*topcs); % projection matrix       
        y=x(:,oChan)*(topcs*b'); % projection 
                
        dx=abs(y-x(:,iChan));   % difference from projection
        dx=dx+eps;              % avoids error on simulated data       
        d=dx/mean(dx(find(w))); % excentricity measure
        d(find(isnan(d)))=0;
        if NSMOOTH>0; 
            d=filtfilt(ones(NSMOOTH,1)/NSMOOTH,1,d);
        end
       
        thresh1=thresh;     % threshold for shared weight, for estimating variance
        thresh2=thresh/THRESH_RATIO;   % threshold for channel-specific weight, to define clean subspace

        w=min(w,d<thresh1);     % w==0 for artifact part
        ww(:,iChan)=d<thresh2;  % weights specific to each channel
        
   end    
    
    prop=mean(w);
    if VERBOSE>0; disp(['proportion artifact free: ', num2str(prop)]); end    
    if iIter==NITER && prop<MINPROP
        thresh=thresh*1.1;
        if VERBOSE>0; disp(['Warning: nt_star increasing threshold to ', num2str(thresh)]); end
        w=ones(size(w));
    else
        iIter=iIter-1;
    end
    
    x=nt_demean(x,w);
    c0=nt_cov(x,[],w);      % recalculate covariance on non-artifact part
    w=ones(size(x,1),1);    % reset
end
ww=nt_growmask(ww,10);


%{
We now know which part contains channel-specific artifacts (w==0 for artifact part), 
and we have an estimate of the covariance matrix of the artifact-free part.
%}

figure(1); clf; plot(mean(ww,2)); 
y=nt_inpaint(x,ww,w);

x=y;


x=nt_demean(x);
x=bsxfun(@times,x,nn);
x=bsxfun(@plus,x,mn);
x=nt_fold(x,nsample);
w=nt_fold(w,nsample);
ww=nt_fold(ww,nsample);
x(:,iNan)=nan;
x(:,iZero)=0;
if VERBOSE>0; disp(['power ratio: ', num2str(nt_wpwr(x)/p0)]); end

function [ww,nOccurrences,iBack]=patternDict(w);
% ww: dictionary of patterns
% nOccurrences: number of times each pattern occurred
% iBack: index to reconstruct input from dictionary
[ww,~,IC,nOccurrences]=nt_unique(w,'rows');
[nOccurrences,iSort]=sort(nOccurrences, 'descend'); % sort by decreasing number
[~,iReverse]=sort(iSort); % 
ww=ww(iSort,:); % same order for patterns, w = ww(iReverse1(IC),:)
iBack=iReverse(IC); % w = ww(iBack,:)

function [w,x]=coalesceWeights(w,x)
% reduce the number of weight patterns based on pruning & interpolation, to
% save time
nchans=size(w,2);
w1=w(1:end-2,:); 
w2=w(2:end-1,:); 
w3=w(3:end,:); 
x1=x(1:end-2,:); 
x2=x(2:end-1,:); 
x3=x(3:end,:); 
% weights that correspond to a single sample 
idx=find(~w1&w2&~w3); % isolated good, set to bad
w2(idx)=0; w(2:end-1,:)=w2;
idx=find(w1&~w2&w3); % isolated bad, interpolate, set to good
w2(idx)=1; w(2:end-1,:)=w2; 
x2(idx)=(x1(idx))+x3(idx)/2; x(2:end-1,:)=x2;
% weight patterns that occur rarely
[wDict,nOccurrences,iBack]=patternDict(w); 
idx=find(nOccurrences<=3);
isolated=iBack(idx);
for k=1:numel(isolated)
    iIsolated=isolated(k);
    if iIsolated-1>=1 && iIsolated+1<=size(w,1)            
        if sum(abs(w(iIsolated-1,:)-w(iIsolated,:)))<sum(abs(w(iIsolated,:)-w(iIsolated+1,:)));
            iFix=iIsolated-1;% preceding is closest
        else
            iFix=iIsolated+1;% following is closest
        end
        idx=find(w(iIsolated,:)&~w(iFix,:));
        w(iIsolated,idx)=0;
        idx=find(~w(iIsolated,:)&w(iFix,:));
        w(iIsolated,idx)=1;
        x(iIsolated,idx)=(x(iIsolated-1,idx)+x(iIsolated+1,idx))/2;
    end
end

    


%% test code
if 0
    %[p,x]=nt_read_data('/data/meg/theoldmanandthesea/eeg/BQ/BQ_aespa_speech_47.mat');
    %[p,x]=nt_read_data('/data/meg/theoldmanandthesea/eeg/MC/MC_aespa_speech_14.mat');
    [p,x]=nt_read_data('/data/meg/theoldmanandthesea/eeg/NM/NM_aespa_speech_12.mat');
    x=x'; x=x(:,1:128); x=x(:,:);
    x=nt_demean(x);
    [x,w]=nt_detrend(x,10); 
    iBad=nt_badChannels(x,[],5000); 
    if ~isempty(iBad)
        disp(['bad channels: ', num2str(iBad)]);
        w(:,iBad)=0;
    end
    
    [xx,w]=nt_star2(x,3,[],w); 
    figure(1); clf; 
    subplot 311; plot(x); title('raw'); subplot 312; plot(xx); title('clean'); subplot 313; plot (x-xx); title('raw-clean');
    ch=1; figure(2); clf; subplot 211; plot ([x(:,ch),xx(:,ch)]); title('single channel'); legend('raw','clean'); subplot 212; plot (x(:,ch)-xx(:,ch)); title('raw-clean')
end   

        

