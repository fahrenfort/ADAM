function [idx,score_vector,score]=nt_split(x,depth,thresh,guard,minstep)
%[idx,score_vector,score]=nt_split(x,depth,thresh,guard,minstep) - split time series into intervals
%
%  idx: index at which to split
%  score: proportion variance reduced
%  score_vector: score as a function of split position
%  
%
%  x: data (time * channels)
%  depth: recursion depth (we find up to 2^depth-1 split points).
%  thresh: if relative variance reduction smaller don't spit [default: 0]
%  guard: number of samples to avoid at each end (counteracts bias towards splitting near ends)
%  minstep: if absolute local step size smaller don't split
%  
%
%  This routine finds the best place to split a time series.
%  
%
% Examples: 
%   nt_split(x);            % find point of largest change
%   nt_split(x.^2);         % largest change of variance
%   nt_split(nt_xprod(x);   % largest change of covariance
%   nt_split(nt_xprod(x,'nodiag'); % same, ignoring variance
%   nt_split(nt_xprod(nt_normrow(x),'nodiag'); % same, each slice normalized
%   nt_split(x,3);          % recurse 3 times (--> 7 split points)
%   nt_split(x,3,10);       % same, but splits are at least 10 points from ends
%   test_nt_split;          % test script

if nargin<2||isempty(depth); depth=1; end
if nargin<3||isempty(thresh); thresh=0; end
if nargin<4||isempty(guard); guard=0; end
if nargin<5; minstep=[]; end

if ndims(x)>2; error('!'); end
if numel(depth)>1; error('!'); end

[m,n]=size(x);
if m<2*guard; % don't try to split if too small
    idx=[]; score_vector=[]; 
    return; 
end

x=nt_demean(x);

%{
For each potential split point we calculate the sum of the per-interval
ssq of first and second interval. This is vectorized using cumsum.
%}

% to minimize memory requirements code is repeated after flipping:
xxx=x;
first_term = cumsum(xxx.^2) - bsxfun(@times, cumsum(xxx).^2,1./(1:m)');
xxx=flipud(x); 
second_term = cumsum(xxx.^2) - bsxfun(@times, cumsum(xxx).^2, 1./(1:m)'); %clear x
score_vector=first_term+second_term(end:-1:1,:);    % sum per-interval ssqs
score_vector=score_vector*diag(1./sum(xxx.^2));     % normalize each dimension
score_vector=mean(score_vector,2);                  % average across dimensions

% find the sweet spot:
[score0,idx]=min(score_vector);

if score0>1-thresh ...                          % improvement too small
        || idx<guard || idx>size(x,1)-guard     % too close to edge
    idx=[]; 
    score_vector=[];
end

if depth>1 && ~isempty(idx)
    [a,sv]=nt_split(x(1:idx,:),depth-1,thresh,guard,minstep);
    [b,sv]=nt_split(x(idx+1:end,:),depth-1,thresh,guard,minstep);
    idx=[a,idx,idx+b];
    idx=unique(idx);
end

% weed out splits with small variance reduction score, or small step size
toWeed=[];
idx=[0, idx, size(x,1)];
for iSplit=1:numel(idx)-2;
    a=nt_demean(x(idx(iSplit)+1:idx(iSplit+1),:));
    b=nt_demean(x(idx(iSplit+1)+1:idx(iSplit+2),:));
    c=nt_demean(x(idx(iSplit)+1:idx(iSplit+2),:));
    MARGIN=10;                              % samples, distance on each side at which to sample step size
    if (sum(a(:).^2)+sum(b(:).^2)) > (1-thresh)*sum(c(:).^2)...
        || (~isempty(minstep) && ...        % local step size is...
        mean(abs(x(idx(iSplit+1)-20,:) - x(idx(iSplit+1)+20,:)),2) < minstep); % ...too small
        toWeed=[toWeed,iSplit];
    end
end
idx(toWeed)=[];

if depth==6
    % refine split positions (sometimes they are off by one or two samples)
    idx2=idx;
    disp(idx)
    for iSplit=1:numel(idx)-2
        ii=nt_split( x(idx(iSplit)+1 : idx(iSplit+2)-2), 1); % -2 to avoid influence of next step
        disp([idx(iSplit) idx(iSplit+2) ii])
        idx2(iSplit+1)=idx(iSplit)+ii;
    end
    idx=idx2;
end

idx=idx(2:end-1);


            

if nargout>2 || nargout==0;
    % score = reduction in variance if mean removed from each segment
    ssq_total=sum( (x-repmat(mean(x),size(x,1),1)).^2 );
    idx2=unique([1,idx,m]); 
    ssq=zeros(1,size(x,2));
    for iSegment=1:numel(idx2)-1
        xx=x(idx2(iSegment):idx2(iSegment+1),:);
        ssq=ssq + sum( (xx - repmat(mean(xx),size(xx,1),1) ).^2 );
    end
    score=(mean(ssq))./mean(ssq_total);
end

%disp(['nt_split_nargout: ', num2str(nargout)])

if nargout==0;
    disp(['split at ', num2str(idx)]);
    disp(['(%: ', num2str(100*idx/m, '  %.01f'), ')'])
    disp(['score: ', num2str(score,  '%.01f')]);
    
    figure(200);
    subplot 211
    plot(score_vector);
    subplot 212
    plot(x); drawnow
    nt_mark(idx);
    if numel(idx)>1; disp(['smallest interval: ', num2str(min(diff(idx)))]); end
    clear idx score score_vector
end

%disp([depth, idx])

% TODO:
% - online "bottom-up" calculation (aggregate rather than split)



