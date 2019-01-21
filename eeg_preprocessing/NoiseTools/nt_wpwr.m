function [y,tweight]=nt_wpwr(x,w)
%[y,tweight]=nt_wpwr(x,w) - weighted power
%
%  y: weighted ssq of x
%  tweight: total weight
%
%  x: data
%  w: weight
%

if nargin<2; w=[]; end

x=nt_unfold(x);
w=nt_unfold(w);

if isempty(w)
    y=sum(x(:).^2);
    tweight=numel(x);
else
    x=nt_vecmult(x,w);
    y=sum(x(:).^2);
    tweight=sum(w(:));
end

if nargout==0; 
    disp(num2str(y));
    clear y tweight
end