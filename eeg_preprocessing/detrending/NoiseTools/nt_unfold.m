function x=nt_unfold(x)
%y=nt_fold(x) - unfold 3D to 2D 
%
%  y: 2D matrix of concatentated data (time * channel)
%
%  x: 3D matrix of (time * channel * trial)
%
% NoiseTools
nt_greetings;


if isempty(x)
    x=[];
else
    [m,n,p]=size(x);
    if p>1;
        x=reshape(permute(x,[1 3 2]), m*p,n);
    else
        x=x;
    end
end