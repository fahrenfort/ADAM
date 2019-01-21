function [iMorton,toMorton]=nt_morton(nrows,ncols)
%[iMorton,toMorton]=nt_morton(nrows,ncols) - indices for Morton scan of image
% 
%  iMorton: matrix of indices in Morton order
%  toMorton: indices to put pixels in Morton order
%
%  nrows,ncols: dimensions of image
%
% NoiseTools
nt_greetings;

if nargin<2; error('!'); end

% morton order for square of side N power of 2
N=2^ceil(log2(max(nrows,ncols)));
if N==2
   iMorton = [1 2; 3 4];
else
   b = nt_morton(N/2,N/2);
   iMorton = [b b+(N/2)^2; b+(N/2)^2*2 b+(N/2)^2*3];
end

% clip to bounds
iMorton=iMorton(1:nrows,1:ncols);

% adjust to avoid skipped values
[~,idx]=sort(iMorton(:));
iMorton(idx)=1:numel(iMorton);

[~,toMorton]=sort(iMorton(:));
