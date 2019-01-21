function sgn=nt_peaksign(x,dim)
%sgn=peaksign(x,dim) - sign of largest extremum
% 
% If dim is present, give array of signs along that dimension

if nargin<2; dim=[]; end

if isempty(dim)
    if max(-x(:))>max(x(:));
        sgn=-1;
    else
        sgn=1;
    end
else
    if dim>4; error('1'); end
    [m,n,o,p]=size(x);
    if dim==1
        for k=1:m
            sgn(k)=nt_peaksign(x(k,:,:,:));
        end
    elseif dim==2
        for k=1:n
            sgn(k)=nt_peaksign(x(:,k,:,:));
        end
    elseif dim==3
        for k=1:o
            sgn(k)=nt_peaksign(x(:,:,k,:));
        end
    else
        for k=1:p
            sgn(k)=nt_peaksign(x(:,:,:,k));
        end
    end
end
        
