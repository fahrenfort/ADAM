function x=nt_fixsign(x)
%y=nt_fixsign(x) - flip signs to maximize inter-component correlation

if ndims(x)<3;
    y=nt_pca(x);
    y=y(:,1)'*x; % correlation with PC
    x=x*sign(diag(y));
else
    [m,n,o]=size(x);
    if 0
        % can't remember the logic of this...
        for k=1:n
            x(:,k,:)=reshape(nt_fixsign(squeeze(x(:,k,:))), [m,1,o]);
        end
    else
        x=nt_fold(nt_fixsign(nt_unfold(x)),m);
    end
end
    
% flip whold matrix so highest peak is positive
if abs(max(x(:)))<abs(min(x(:)))
    x=-x;
end