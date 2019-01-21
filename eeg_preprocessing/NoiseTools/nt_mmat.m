function y=nt_mmat(x,m)
%y=nt_mmat(x,m) -  matrix multiplication (with convolution)
%
%  y: result
% 
%  x: input data (2D or more)
%  m: matrix to apply (2D: right multiply, 3D: same with convolution)
%
%  If m is 3D, the last index (k) is lag.

if nargin<2; error('!'); end

if iscell(x)
    for iCell=1:length(x)
        y{iCell}=nt_mmat(x{iCell},m);
    end
    return;
end

if ndims(x)>3
    % concatenate the last dimensions, process, then de-concatenate
    sz=size(x);
    x=reshape(x,[sz(1),sz(2),prod(sz(3:end))]);
    x=nt_mmat(x,m);
    x=reshape(x,[size(x,1),sz(2),sz(3:end)]);

else

    if ndims(m)==2;
        % no convolution
        y=nt_mmat0(x,m); 
    
    else
        
% does anyone use this ?????

        [nRows,nCols,nLags]=size(m);
        [nSamples,nChans,nTrials]=size(x);
        if nChans~=nRows; 
            error('ncols(x) ~= nrows(m)');
        end
        % convolution: for each k, multiply x by m(:,:,k) and add with
        % shift of (k-1)
        y=zeros(nSamples+nLags-1,nCols,nTrials);
        for iLag=1:nLags
            y(iLag:iLag+nSamples-1,:,:) = y(iLag:iLag+nSamples-1,:,:) + nt_mmat0(x,m(:,:,iLag));
        end
        
    end
end    

function x=nt_mmat0(x,m)
x=nt_fold(nt_unfold(x)*m,size(x,1));
