function z=nt_multishift(x,shifts)
%z=nt_multishift(x,shifts,amplitudes) - apply multiple shifts to matrix
%
%   y: result
%
%   x: matrix to shift
%   shifts: array of shifts (must be nonnegative)
% 
% X is shifted column by column (all shifts of 1st column, then all
% shifts of second column, etc).
% 
% X may be 1D, 2D or 3D. See also convmtx.
%
% NoiseTools

if nargin<2; error('!'); end

if iscell(x)
    for iCell=1:length(x);
        z{iCell}=nt_multishift(x{iCell},shifts);
    end
    return;
end

if size(x,1)<max(shifts); error('shifts should be no larger than nrows'); end
if min(shifts)<0; error('shifts should be nonnegative'); end
shifts=shifts(:)';
nshifts=numel(shifts);
if nshifts==1 && shifts(1)==0; 
    z=x;
    return
end

% array of shift indices
N=size(x,1)-max(shifts); 
shiftarray=nt_vecadd(nt_vecmult(ones(N,nshifts),shifts),(1:N)');
[m,n,o]=size(x);
z=zeros(N,n*nshifts,o);

for k=1:o
    for j=0:n-1
        y=x(:,j+1,k);
        z(:,j*nshifts+1: j*nshifts+nshifts,k)=y(shiftarray);
    end
end

