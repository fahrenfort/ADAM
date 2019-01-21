 function [y,yy]=nt_inpaint(x,w)
%function y=nt_inpaint(x,w) - weighted interpolation based on correlation structure
%
%  y: interpolated data (only over w==0, a la STAR)
%  yy: interpolated data (all times, a la SNS)
%
%  x: data matrix
%  w: weight vector or matrix
%
%
% Noisetools.
nt_greetings;

PCA_THRESH=10^-15;
NSMOOTH=0;
thresh=0.5;
nClosest=min(20,size(x,2)-1);
%nClosest=40;

if nargin<1; error('!'); end
if nargin<2||isempty(w); w=ones(size(x)); end
if ndims(x)>2; error('!'); end
if ~all(size(x)==size(w)); error('!'); end
[nsamples,nchan]=size(x);

%{ 
We have multichannel data x and a multichannel weighting function w (0 or
1). There are many configurations of valid/invalid channels to consider.
List them.
%}

[ww,nOccurrences,iBack]=patternDict(w); % all patterns of good/bad channels
nPatterns=size(ww,1);
disp(size(ww));

%{ 
Now we have a list of all the different weight patterns: ww. The
vector iBack indicates which data samples fit each pattern: w = ww(iBack,:).
%}
     
%{
Find which channels are 'neighbors' in terms of covariance.
%}

% weighted covariance matrix to determine which channels are close
x0=x;
[x,save_mean]=nt_demean(x,w); 
[x,save_amp]=nt_normcol(x,w);
xx=x.*w;
c=(xx'*xx) ./ (w'*w); clear xx;
c=abs(c); 
sims=c+10*eye(size(c)); % make sure self always scores highest so we can skip it

y=x;  

%{
We now have a matrix indicating proximity between channels. 
%}

%{
For each channel, we calculate the projection matrix on the the subspace spanned 
by other *valid* channels.  There are as many projection matrices as patterns 
of valid/invalid channels.  Each projection matrix is estimated on data samples for
which iChan is valid, and can be used to reconstruct data samples for which it is 
invalid.
%}

for iChan=1:nchan

    %{ 
    We want to avoid having to consider all patterns of valid/unvalid 
    other channels. For that we'll group patterns. 
    First we order the other channels by decreasing similarity, putting
    invalid samples last. This needs to be done for each pattern.
    %}

    sim=sims(iChan,:);              % list of similarities with all other channels
    sim=repmat(sim,nPatterns,1);    % replicate so each pattern has own list
    sim((~ww))=0;                   % for each list, give bad channels a low score
    [~,closest]=sort(abs(sim),2,'descend');     % sort each list by decreasing similarity
    for k=1:size(closest,1);
        closest(k,find(sim(k,closest(k,:))==0))=0;     % mark bad channels as 0
    end
    for k=1:size(closest,1);
        if closest(k,1)==iChan; 
            closest(k,1:size(closest,2)-1) = closest(k,2:end); % skip first entry of list if same as iChan
        else
            % iChan was bad so not first
        end
    end
    closest=closest(:,1:end-1); % last not valid if first skipped
    
    %{ 
    We now have, for each pattern, a list of channels closest to iChan. 
    There are a lot of different patterns, so we merge those for which the nClosest 
    channels are the same.
    %}
    
    % group patterns for which the nClosest most similar channels are the same
    [C,IA,IC]=unique(closest(:,1:nClosest),'rows');
    iBack2=IC(iBack);       % maps each pattern to the data that fit it
    
    %{
    We now have a smaller array C of reduced patterns. The
    vector iBack2 indicates which data samples correspond to each pattern.
    %}
    
    %{
    For some patterns, iChan is valid throughout.  We can skip these.
    For others, only a few samples are invalid. To save time w can skip these
    too and fix them later using serial interpolation.
    %}
    
    toFix=[];
    NSKIP=2;
    www=ones(size(x,1),1);
    for iPattern=1:size(C,1)
        %if any(~w(find(iBack2==iPattern),iChan)); toFix=[toFix,iPattern]; end
        mySamples=find(iBack2==iPattern); 
        mySamples=mySamples(find(~w(mySamples,iChan)));
        if numel(mySamples)<=NSKIP
            www(mySamples)=0;
        else
            toFix=[toFix,iPattern];
        end
    end
    C=C(toFix,:);
           
    for iPattern=1:size(C,1)

        %{ 
        Estimate matrix to project iChan on the other channels listed in this
        pattern. 
        %}

        oChan=C(iPattern,:);
        oChan(find(oChan==0))=[]; % exclude bad channels

        % samples corresponding to this pattern
        mySamples=find(iBack2==toFix(iPattern)); 
        mySamples=mySamples(find(~w(mySamples,iChan)));
%        mySamples=find(iBack2==iPattern); 
                
        % select data for which iChan *and* oChan are valid
        iBothValid=all(w(:,[iChan,oChan]),2);        
        xxx=x(iBothValid, [iChan,oChan]);  
        
        %if size(xxx,1)<8000; disp([iChan, iPattern]); disp(size(xxx,1)); end


        %%% --> we should be able to avoid this situation
        if isempty(xxx); 
            disp([iChan, iPattern]); disp('empty'); 
            continue; % we won't estimate or fix anything
        end
        
        % calculate covariance matrix
        mn=mean(xxx,1);
        xxx=nt_demean(xxx); % remove mean first
        ccc=xxx'*xxx;

        % PCA other channels to remove weak dimensions
        [topcs,eigenvalues]=nt_pcarot(ccc(2:end,2:end));
        idx=find(eigenvalues/max(eigenvalues) > PCA_THRESH); % discard weak dims
        topcs=topcs(:,idx);

        % projection matrix
        b=ccc(1,2:end)*topcs / (topcs'*ccc(2:end,2:end)*topcs); 
        
        %{
        We now have the projection matrix to project channel iChan on channels oChan,
        applicable to samples corresponding to this pattern.  We can use it
        to fix samples for which iChan is invalid.
        %}
        
        y(mySamples,iChan) = ...
            (x(mySamples,oChan) - repmat(mn(2:end),numel(mySamples),1))... % first remove mean of other chans...
            *(topcs*b') ...
            + mn(1); % ... then restore mean of this channel
        
    end
    
    %{
    Now we fix the isolated samples that we skipped using serial interpolation.
    %}
    MAXGAPSIZE=100;
    y(:,iChan)=fillgap(y(:,iChan),www,MAXGAPSIZE);
    
end


y=bsxfun(@times,y,save_amp); % restore the initial ampitude
y=bsxfun(@plus,y,save_mean); % restore the initial mean

yy=y; 
y(w~=0) = x0(w~=0); % don't touch valid parts

if ~nargout
    % plot, don't return values
    disp(nt_wpwr(y)/nt_wpwr(x));
    figure(11); clf;
    subplot 311; plot(x0); title('raw');
    subplot 312; plot(y); title('projected on other channels')
    subplot 313; plot(w); ylim([-.1 1.1]); title('weight');
    clear w y
end

function [y,w]=fillgap(x,w,maxGapSize)
%y=fillgap(x,w,maxGapSize) - fill gaps using simple interpolation
%
%  y: interpolated data
% 
%  x: data to interpolate
%  w: weighting function
%  maxGapSize: largest expected gap size
if nargin<2; error('!'); end
if nargin<3||isempty(maxGapSize); maxGapSize=1; end
if size(x,2)>1; error('!'); end
if size(x) ~= size(w); error('!'); end
y=x;
if all(w); return; end
% simple case size one
iToFix=1+find(~w(2:end-1)&w(1:end-2)&w(3:end)); 
y(iToFix)=(y(iToFix-1)+y(iToFix+1))/2;
w(iToFix)=1; 
% general case size > 1
iStart=find(w(1:end-2) & ~w(2:end-1));  % position preceding gap
iStop=find(~w(1:end-1) & w(2:end));     % last position in gap
if isempty(iStart)||isempty(iStop); return; end
if iStop(1)<iStart(1);
    iStop=iStop(2:end);                 % ignore gap at beginning
end
iStart=iStart(1:numel(iStop));          % ignore gap at end
for gapSize=2:maxGapSize
    idx=find(iStop-iStart==gapSize);
    for k=1:gapSize
        % interpolate between samples on either side of gap
        y(iStart(idx)+k) = ( y(iStart(idx)) * (gapSize-k+1) + y(iStart(idx)+gapSize+1) * k ) / (gapSize+1);
        w(iStart(idx)+k) = 1;
    end
end


% create a dictionary of weight patterns
function [ww,nOccurrences,iBack]=patternDict(w)
% ww: dictionary of patterns
% nOccurrences: number of times each pattern occurred
% iBack: index to reconstruct input from dictionary
[ww,~,IC,nOccurrences]=nt_unique(w,'rows');
[nOccurrences,iSort]=sort(nOccurrences, 'descend'); % sort by decreasing number
[~,iReverse]=sort(iSort); % 
ww=ww(iSort,:); % same order for patterns, w = ww(iReverse1(IC),:)
iBack=iReverse(IC); % w = ww(iBack,:)

%%% TEST %%%
if 0
    x0=sin(2*pi*(1:10000)'*(1:5)/10000);
    x=x0*randn(5,10)+1;
    x(1:4000,1)=x(1:4000,1)+0.3*randn(size(x(1:4000,1)));
    w=ones(size(x));  
    w(1:4000,1)=0; w(4001:6000,3)=0; %w(6001:7000,4)=0; w(7001:8000,5)=0; w(8001:9000,6)=0; w(9001:10000,7)=0;
    %w=randn(size(w))>0;
    %b=nt_regw(x(:,1),x(:,2:end),w(:,1)); y=x(:,2:end)*b;
    y=nt_inpaint(x,w);
    figure(2); clf
    subplot 211; plot([x(:,1),y(:,1)]); legend('raw','clean'); title('one channel'); subplot 212; plot(x(:,1)-y(:,1)); title('raw-clean');
end
if 0
    x0=sin(2*pi*(1:10000)'*(1:3)/10000);
    x=x0*nt_normcol(randn(3,5));
    w=ones(size(x));
    x(1:1000,1)=100; w(1:1000,1)=0;
    x(2001:3000,1)=100; w(2001:3000,1)=0;
    x(1:2000,2)=100; w(1:2000,2)=0;
    x=x+0.1*randn(size(x));
    [y,yy]=nt_inpaint(x,w);
    figure(1); clf
    subplot 311; plot(x); title('raw');  subplot 312; plot(y); title('clean'); subplot 313; plot(x-y); title('raw-clean');
end
if 0
    N=3;
    nchans=50;
    x=zeros(1100,N);
    for k=1:N
        x(:,k)=sin(2*pi*k*(1:1100)'/1100);
    end
    x=x*randn(N,nchans);
    x=nt_normcol(x) + 0*randn(size(x));
    xx=x;
    w=ones(size(x));
    for k=1:nchans
        xx(k*20+(1:10),k)=100;
        w(k*20+(1:10),k)=0;
    end
    [y,yy]=nt_inpaint(xx,w);
    figure(1); clf;
    subplot 411; plot(x); title('orig'); subplot 412; plot(xx); title('w glitches'); 
    subplot 413; plot(y); title('clean'); subplot 414; plot(x-y); title('diff w orig');
end
if 0
    N=10;
    nchans=20;
    nsamples=1100;
    x=zeros(nsamples,N);
    for k=1:N
        x(:,k)=sin(2*pi*k*(1:nsamples)'/nsamples);
    end
    x=x*randn(N,nchans);
%    x=x+1*randn(size(x)); % add noise
    xx=x;
    w=ones(size(x));
    for k=1:nchans
        xx(500+k*20+(1:40),k)=100;
        w(500+k*20+(1:40),k)=0;
    end
    [y]=nt_inpaint(xx,w);
    figure(1); clf;
    subplot 411; plot(x); title('original'); 
    subplot 412; plot(xx); title ('with glitches');
    subplot 413; plot (y); title ('fixed'); 
    subplot 414; plot(x-y); title ('error');
end
if 0
    [p,x]=nt_read_data('/data/meg/theoldmanandthesea/eeg/mc/MC_aespa_speech_45.mat'); x=x'; x=x(:,1:128); x=x(0+(1:10000),:);
    %[p,x]=nt_read_data('/data/meg/arzounian/ADC_DA_140521_p20/ADC_DA_140521_p20_01_calib'); x=x'; x=x(1:10000,:);
    
    x=nt_demean(x);
    [x,w]=nt_detrend(x,10);   
    profile on; y=nt_inpaint(x,w); profile report;
    figure(1); clf
    subplot 311; plot(x); title('raw');  subplot 312; plot(y); title('clean'); subplot 313; plot(x-y); title('raw-clean');
    figure(2); clf
    ch=35;subplot 311; plot([x(:,ch),y(:,ch)]); subplot 312; plot(x(:,ch)-y(:,ch)); subplot 313; plot(w(:,ch), '.-');
end

    