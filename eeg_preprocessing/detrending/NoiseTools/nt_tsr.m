function [y,idx,w]=nt_tsr(x,ref,shifts,wx,wref,keep,thresh)
%[y,idx,w]=nt_tsr(x,ref,shifts,wx,wref,keep,thresh) - time-shift regression (TSPCA)
%
%  y: denoised data 
%  idx: x(idx) is aligned with y
%  w: weights applied by tsr
% 
%  x: data to denoise (time * channels * trials)
%  ref: reference (time * channels * trials)
%  shifts: array of shifts to apply to ref (default: [0])
%  wx: weights to apply to x (time * 1 * trials);
%  wref: weights to apply to ref (time * 1 * trials);
%  keep: number of shifted-ref PCs to retain (default: all)
%  thresh: ignore shifted-ref PCs smaller than thresh (default: 10.^-12)
%
% NoiseTools
nt_greetings;


% Copyright 2007, 2008 Alain de Cheveigne

% See: 
% de Cheveign\'e, A. and Simon, J. Z. (2007). "Denoising based on
% Time-Shift PCA." Journal of Neuroscience Methods 165: 297-305.
%
% The basic idea is to project the signal X on a basis formed by the
% orthogonalized time-shifted REF, and remove the projection. Supposing REF
% gives a good observation of the noise that contaminates X, the noise is
% removed. By allowing time shifts, the algorithm finds the optimal FIR filter 
% to apply to REF so as to compensate for any convolutional mismatch
% between X and REF.
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
%

if nargin<2; error('too few arguments'); end
if nargin<3 || isempty(shifts); shifts=0; end
if nargin<4; wx=[]; end
if nargin<5; wref=[]; end
if nargin<6 || isempty(keep); keep=[]; end
if nargin<7 || isempty(thresh); thresh=10.^-20; end

% check argument values for sanity
if size(x,1)~=size(ref,1); error('X and REF should have same nrows'); end
if size(x,3)~=size(ref,3); error('X and REF should have same npages'); end
if ~isempty(wx) && size(x,3)~=size(wx,3); error('X and WX should have same npages'); end
if ~isempty(wx) && size(x,1)~=size(wx,1); error('X and WX should have same nrows'); end
if ~isempty(wref) && size(ref,1)~=size(wref,1); error('REF and WREF should have same nrows'); end
if ~isempty(wref) && size(ref,3)~=size(wref,3); error('REF and WREF should have same npages'); end
if max(shifts)-min(0,min(shifts)) >= size(x,1); error('X has too few samples to support SHIFTS'); end
if ~isempty(wx) && size(wx,2)~=1; error('wx should have ncols=1'); end
if ~isempty(wref) && size(wref,2)~=1; error('wref should have ncols=1'); end
if ~isempty(wx) && sum(wx(:))==0; error('weights on x are all zero!'); end
if ~isempty(wref) && sum(wref(:))==0; error('weights on ref are all zero!'); end

if numel(shifts)>1000; error(['numel(shifts)=',num2str(numel(shifts)), ' (if OK comment out this line)']); end


% We need to adjust x and ref to ensure that shifts are non-negative.  
% If some values of shifts are negative, we increment shifts and truncate x.

% adjust x to make shifts non-negative
offset1=max(0,-min(shifts));
idx=1+offset1:size(x,1);
x=x(idx,:,:);                             % truncate x
if ~isempty(wx); wx=wx(idx,:,:); end
shifts=shifts+offset1;                    % shifts are now positive

% adjust size of x
offset2=max(0,max(shifts)); 
idx=1: size(ref,1)-offset2; 
x=x(idx,:,:);                           % part of x that overlaps with time-shifted refs
if ~isempty(wx); wx=wx(idx,:,:); end

[mx,nx,ox]=size(x);
[mref,nref,oref]=size(ref);

% consolidate weights into single weight matrix
w=zeros([mx,1,oref]);
if isempty(wx) && isempty(wref)
    w(1:mx,:,:)=1;
elseif isempty(wref);
    w(:,:,:)=wx(:,:,:);
elseif isempty(wx)
    for k=1:ox
        wr=wref(:,:,k);
        wr=ts_multishift(wr,shifts);
        wr=min(wr,[],2);
        w(:,:,k)=wr;
    end;
else
    for k=1:ox
        wr=wref(:,:,k);
        wr=nt_multishift(wr,shifts);
        wr=min(wr,[],2);
        wr=min(wr,wx(1:size(wr,1),:,k));
        w(:,:,k)=wr;
    end
end
wx=w;
wref=zeros(mref,1,oref);
wref(idx,:,:)=w;

% remove weighted means
x0=x;
x=nt_demean(x,wx);
mn1=x-x0;
ref=nt_demean(ref,wref);

% equalize power of ref chans, then equalize power of ref PCs
ref=nt_normcol(ref,wref);
ref=nt_pca(ref,0,[],10^-6);
ref=nt_normcol(ref,wref);

% covariances and cross covariance with time-shifted refs
[cref,twcref]=nt_cov(ref,shifts,wref);
[cxref,twcxref]=nt_xcov(x,ref,shifts,wx);

% regression matrix of x on time-shifted refs
r=nt_regcov(cxref/twcxref,cref/twcref,keep,thresh);

%r=r*0.765;

% TSPCA: clean x by removing regression on time-shifted refs
y=zeros(mx,nx,ox);
for k=1:ox
    z=nt_multishift(ref(:,:,k),shifts)*r;
    %plot([x(1:size(z,1),1,k), z(:,1)]); pause
    y(:,:,k)=x(1:size(z,1),:,k)-z;
end
y0=y;
y=nt_demean(y,wx);    % multishift(ref) is not necessarily 0 mean
mn2=y-y0;

%idx=1+offset1:n0-offset2;
idx=1+offset1:size(y,1)+offset1;
mn=mn1+mn2;
w=wref;

%y=vecadd(y,mn);
