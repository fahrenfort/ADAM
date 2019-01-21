function [x,w,ww]=nt_star(x,thresh,closest,depth)
% [y,w,ww]=nt_star(x,thresh,closest,depth) - sensor noise suppression
%
%  y: denoised data 
%  w: 0 for parts that needed fixing, 1 elsewhere (time*1)
%  ww: 0 for parts that needed fixing, 1 elsewhere (time*chans)
%
%  x: data to denoise (time*chans or time*chans*trials)
%  thresh: threshold for excentricity measure (default:1);
%  closest: indices of channels that are closest to each channel (default: all)
%  depth: maximum number of channels to fix at each sample (default 1)
% 
% See also: nt_sns, nt_proximity

nt_greetings;

PCA_THRESH=10^-15;  % threshold for discarding weak PCs
NSMOOTH=10;         % samples, smoothing to apply to excentricity
MINPROP=0.3;        % minimum proportion of artifact free at first iteration
NITER=3;            % iterations to refine c0
VERBOSE=1;          % set to 0 to shut up

if nargin<1; error; end
if nargin<2 || isempty(thresh); thresh=1; end
if nargin<3; closest=[]; end
if ~isempty(closest)&&size(closest,1)~=size(x,2);
    error('closest array should have as many rows as channels of x'); 
end
if nargin<4 || isempty(depth); depth=1; end

if nargout==0 % plot, don't return result
    [y,w,ww]=nt_star(x,thresh,closest,depth);
    disp([mean(w(:)), mean(ww(:))])
    figure(1); clf;
    subplot 311; plot(x);
    subplot 312; plot(y);
    subplot 313; plot(w, '.-'); ylim([0 1.1]);
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

% NaN and zero channels are set to rand, which effectively excludes them
iNan=find(all(isnan(x)));
iZero=find(all(x==0));
x(:,iNan)=randn(size(x,1),numel(iNan));
x(:,iZero)=randn(size(x,1),numel(iZero));

x=nt_demean(x);
c0=nt_cov(x); % initial covariance estimate

%{
Find time intervals where at least one channel is excentric --> w==0.
%}

iIter=NITER;
while iIter>0
    
    
    w=ones(size(x,1),1);
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
        if NSMOOTH>0; 
            d=filtfilt(ones(NSMOOTH,1)/NSMOOTH,1,d);
        end
        
        d=d/thresh;
        w=min(w,(d<1)); % w==0 for artifact part
        
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
    c0=nt_cov(x,[],w); % restrict covariance estimate to non-artifact part
end

%{
We now know which part contains channel-specific artifacts (w==0 for artifact part), 
and we have an estimate of the covariance matrix of the artifact-free part.
%}

%{
Find which channel is most excentric at each time point.
Here we use an excentricity measure based on the absolute value of the signal,
rather than the difference between signal and projection.
%}

xx=abs(x);
xx=bsxfun(@times,xx, 1 ./ sqrt(mean(xx(find(w),:).^2))); % divide by std over non-artifact part
if NSMOOTH>0; 
    xx=filtfilt(ones(NSMOOTH,1)/NSMOOTH,1,xx);
end
[~,rank]=sort(xx','descend'); 
rank=rank';
rank(find(w),:)=0;      % exclude parts that were not judged excentric

depth=min(depth,nchan-1);
ww=ones(size(x));
for iDepth=1:depth

    %{
    Fix each channel by projecting on other channels.
    %}
    
    iFixed=nchan;
    nFixed=0;
    for iChan=1:nchan

        bad_samples=find(iChan==rank(:,iDepth)); % samples where this channel is the most excentric
        if iDepth ~=1; 
            bad_samples(find(xx(bad_samples,iChan)<thresh)) =[]; % exclude if not very bad            
        end
        
        nFixed=nFixed+numel(bad_samples);
        if isempty(bad_samples); 
            iFixed=iFixed-1;
            continue;
        end
        ww(bad_samples,iChan)=0;

        % other channels
        if ~isempty(closest); 
            oChan=closest(iChan,:);
        else
            oChan=setdiff(1:nchan,iChan);
        end
        oChan(oChan>nchan)=[]; % in case closest includes channels not in data

        % PCA other channels to remove weak dimensions
        [topcs,eigenvalues]=nt_pcarot(c0(oChan,oChan)); % PCA
        idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
        topcs=topcs(:,idx);

        % project this channel on other channels
        b=c0(iChan,oChan)*topcs/(topcs'*c0(oChan,oChan)*topcs); % projection matrix 
        y=x(bad_samples,oChan)*(topcs*b'); % projection 

        x(bad_samples,iChan)=y(:); % fix

    end
    
    if VERBOSE>0; 
        disp(['depth: ', num2str(iDepth), ', n fixed channels: ',num2str(iFixed),...
            ', n fixed samples: ', num2str(nFixed), ', ratio: ',num2str(nt_wpwr(x)/p00)]);
    end
end

x=nt_demean(x);
x=bsxfun(@times,x,nn);
x=bsxfun(@plus,x,mn);

x=nt_fold(x,nsample);
w=nt_fold(w,nsample);
ww=nt_fold(ww,nsample);




x(:,iNan)=nan;
x(:,iZero)=0;

if VERBOSE>0; disp(['power ratio: ', num2str(nt_wpwr(x)/p0)]); end