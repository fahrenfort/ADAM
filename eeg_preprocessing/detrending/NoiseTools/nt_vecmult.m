function x=nt_vecmult(x,v)
%y=nt_vecmult(x,v) - multiply all rows or columns of matrix by vector
%
% See vecadd, bsxfun
%
% NoiseTools

% check once and for all to save time
persistent bsxfun_exists;
if isempty(bsxfun_exists); bsxfun_exists=(exist('bsxfun')==5); end

[m,n,o]=size(x);
x=nt_unfold(x);
v=nt_unfold(v);

[mm,nn]=size(x);
[mv,nv]=size(v);
if mv==mm
    % same number of rows, v should be column vector (or same size as x)
    if nv==nn
        x=x.*v;
    elseif nv==1
        if bsxfun_exists;
            x=bsxfun(@times,x,v);
            %y=vecop_core(x, v, 1, 2);  % 2 is the opcode of multiplication in vecop_core
        else
            x=x .* repmat(v,1,nn);
        end
    else
        error('V should be row vector'); 
    end

elseif nv==nn
    % same number of columns, v should be row vector (or same size as x)
    if mv==mm
        x=x.*v;
    elseif mv==1
        if bsxfun_exists;
            x=bsxfun(@times,x,v);
            %y=vecop_core(x, v, 2, 2);  % 2 is the opcode of multiplication in vecop_core
        else
            x=x .* repmat(v,mm,1);
        end
    else
        error('V should be column vector'); 
    end    

else
    error('V and X should have same number of rows or columns'); 
end

x=nt_fold(x,m);

