function z=nt_multismooth(x,smooth,alignment)
%z=nt_multismooth(x,smooth,alignment) - apply multiple smoothing kernels
%
%   y: result
%
%   x: matrix to smooth
%   smooth: vector of smoothing kernel sizes
%   alignment: -1: left [default], 0: center, 1:right
% 
% X is smoothed column by column (all smoothed versions of 1st column, then all
% of second column, etc).
% 
% X may be 1D, 2D or 3D. See also nt_multishift.
%
% NoiseTools
nt_greetings;

if nargin<3 || isempty(alignment); alignment=-1; end
if nargin<2; error('!'); end
if min(smooth)<1; error('smooth must be positive'); end

if iscell(x)
    for iCell=1:length(x);
        z{iCell}=nt_multismooth(x{iCell},smooth,alignment);
    end
    return
end

if size(x,1)<max(smooth); error('smoothing kernel size should be no larger than nrows'); end
if min(smooth)<0; error('smoothing kernel size should be nonnegative'); end
smooth=smooth(:)';
nsmooth=numel(smooth);

% array of shift indices
[m,n,o]=size(x);
z=zeros(m,n*nsmooth,o);

for iPage=1:o
    zz=zeros(m,n,nsmooth);
    for iSmooth=1:nsmooth
        if alignment==-1; nodelayflag=0; elseif alignment==0; nodelayflag=1; else; error('not implemented'); end
        zz(:,:,iSmooth)=nt_smooth(x(:,:,iPage),smooth(iSmooth),[],nodelayflag);
    end
    zz=permute(zz,[1,3,2]); 
    zz=reshape(zz,m,n*nsmooth);
    z(:,:,iPage)=zz;
end


