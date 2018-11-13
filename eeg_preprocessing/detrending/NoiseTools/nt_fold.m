function x=fold(x,epochsize)
%y=fold(x,epochsize) - fold 2D to 3D 
%
%  y: 3D matrix of (time * channel * trial)
%
%  x: 2D matrix of concatentated data (time * channel)
%  epochsize: number of samples in each trial
%
% NoiseTools

nt_greetings;

if isempty(x); 
    x=[]; 
else
    if size(x,1)/epochsize>1
        x=permute(reshape(x,[epochsize, size(x,1)/epochsize, size(x,2)]), [1 3 2]);
    else
        x=x;
    end
end