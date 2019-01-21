function i=nt_idx(x,scale,i)
%i=nt_idx(x,scale,i) - index a data matrix
%
%  i: index structure
%  
%  x: data 
%  scale: scale 
%  i: already populated structure
%
% The 'i' argument can be used to control which statistics are calculated, e.g.
% i.min=[]; i.max=[]; i=nt_idx(x,n,i); % creates min and max indexes
%
% NoiseTools

if nargin<2; n=100; end
if nargin<3; 
    % default fields
    i.min=[];
    i.max=[];
    i.mean=[];
    i.ssq=[];
    i.card=[];
end

if ndims(x)>2; error('x should be 2D'); end
[m,n]=size(x);
if m<scale; warning('nrows < scale'); end

nchunks=floor(m/scale);
nextra=m-scale*nchunks;

% reshape to calculate stats
extra=x(scale*nchunks+1:end,:);     % extra chunk
x=x(1:scale*nchunks,:);             % main chunks
x=reshape(x,[scale,nchunks,n]);     % reshape to 3D 

% cardinality
if isfield(i,'card')
    a=ones(nchunks,1)*scale;
    if nextra>0; a=[a;nextra]; end
    i.card=[i.card;a];
end

% min
if isfield(i,'min');
    a=reshape(min(x,[],1),[nchunks,n]);
    if nextra>0; a=[a;min(extra,[],1)]; end
    i.min=[i.min;a];
end

% max
if isfield(i,'max');
    a=reshape(max(x,[],1), [nchunks,n]);
    if nextra>0; a=[a;max(extra,[],1)]; end
    i.max=[i.max;a];
end

% sum
if isfield(i,'mean');
    a=reshape(mean(x,1),[nchunks,n]);
    if nextra>0; a=[a;mean(extra,1)]; end
    i.mean=[i.mean;a];
end

% ssq
if isfield(i,'ssq');
    a=reshape(sum(x.^2,1),[nchunks,n]);
    if nextra>0; a=[a;sum(extra.^2,1)]; end
    i.ssq=[i.ssq;a];
end
