 function y=nt_sns(x,nneighbors,skip,w)
% y=nt_sns(x,nneigbors,skip,w) - sensor noise suppression
%
%   y: denoised matrix
%
%   x: matrix  to denoise
%   nneighbors: number of channels to use in projection
%   skip: number of closest neighbors to skip (default: 0)
%   w : weights (default: all ones)
%
% If w=='auto', the weights are calculated automatically.
%
%  Mean of data is NOT removed.
%
% Cite: 
% de Cheveign\'e, A. and Simon, J. Z. (2007). "Sensor Noise Suppression." 
% Journal of Neuroscience Methods, 168: 195-202.
%
% NoiseTools
nt_greetings;


% Copyright 2007, 2008 Alain de Cheveigne


%
% The basic idea is to project each channel of X on a basis formed by the
% orthogonalized set of other channels. Supposing (a) that sensor noise is
% uncorrelated across sensors, and (b) genuine signal is correlated, sensor
% noise is removed and genuine signal preserved. 
% 
% Implementation issues:
% - Data are often available as an array of epochs. This implementation
% caters for 3D data (time * channnels * trials);
% - It is important to deemphasize high amplitude artifacts and glitches
% so that they do not dominate the solution.  This implementation uses
% weighted covariance and means.
% - Processing assumes zero-means data. Means are calculated with weights.
% - The implementation tries to be efficent and minimize memory requirements
% so as to handle large data sets.
%
% Larger data sets (disk based) could be handled by performing mean and
% covariance calculations block-by-block, in several passes.


if nargin<4; w=[]; end
if nargin<3 || isempty(skip); skip=0; end
if nargin<2 || isempty(nneighbors); error('need to specify nneighbors'); end
if ~isempty(w) && sum(w(:))==0; error('weights are all zero!'); end

if ~isempty(find(isnan(x))); error('x contains NANs'); end
if numel(nneighbors)>1 || numel(skip)>1; error('nneighbors & skip should be scalars');  end

[m,n,o]=size(x);
x=nt_unfold(x);

%[x,mn0]=demean(x);  % remove mean
[c,nc]=nt_cov(x);    % raw covariance

TOOBIG=10;
if strcmp(w,'auto')
    y=nt_sns(nt_demean(x),nneighbors,skip);
    d=(y-x).^2;
    d=nt_vecmult(nt_unfold(d), 1./mean( [mean(nt_unfold(x.^2)); mean(nt_unfold(y.^2))] ));
    w=d<TOOBIG;
    w=min(w,[],2);
    w(find(isnan(w)))=0;
end


% sns matrix
if ~isempty(w);
    w=nt_unfold(w);
    %[x,mn1]=demean(x,w);
    [wc,nwc]=nt_cov(nt_demean(x,w),[],w);                     % weighted covariance
    r=nt_sns0(c,nneighbors,skip,wc);
else
    mn1=0;
    w=ones(n,o);
    r=nt_sns0(c,nneighbors,skip,c);
end


% apply to data
y=x*r;

y=nt_fold(y,m);

%mn=mn0;%+mn1;
