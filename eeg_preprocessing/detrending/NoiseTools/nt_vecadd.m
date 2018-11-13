function x=nt_vecadd(x,v)
%y=nt_vecadd(x,v) - add vector to all rows or columns of matrix 
%
% See vecmult, bsxfun
%
% NoiseTools

% check once and for all to save time
persistent bsxfun_exists;
if isempty(bsxfun_exists); 
    bsxfun_exists=(exist('bsxfun')==5); 
    if ~bsxfun_exists; 
        warning('bsxfun not found.  Using repmat');
    end
end


[m,n,o]=size(x);
x=nt_unfold(x);
v=nt_unfold(v);

[mm,nn]=size(x);
if numel(v)==1;
    x=x+v;
elseif size(v,1)==1
    if size(v,2)~=nn; error('V should have same number of columns as X'); end
    if bsxfun_exists;
        x=bsxfun(@plus,x,v);
%        vecop_core(x, v, 2, 1);  % 1 is the opcode of addition in vecop_core
    else
        x=x + repmat(v,mm,1);
    end
elseif size(v,2)==1
    if size(v,1)~=mm; error('V should have same number of rows as X'); end
    if bsxfun_exists;
        x=bsxfun(@plus,x,v);
 %       y=vecop_core(x, v, 1, 1);  % 1 is the opcode of addition in vecop_core
    else
        x=x + repmat(v,1,nn);
    end
end

x=nt_fold(x,m);