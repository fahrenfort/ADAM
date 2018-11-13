function b=nt_lower_to_full(a,ind,n)
% b=nt_lower_to_full(a,ind) - transform lower diagonal to full  covariance
%
%  b: full matrix
%
%  a: matrix of lower diagonal terms
%  ind: indices of lower diagonal terms
%  n: covariance matrix is n*n
%
% Typically used to transform output of nt_xprod into a series of
% covariance matrices.  
%

if nargin<2; error('!'); end
if nargin<3; 
    % estimate based on on size of a (this fails if diagonal is not
    % present)
    n=floor(sqrt(2*size(a,2)));
end

b=zeros(size(a,1),n,n);

[I,J]=ind2sub(n,ind);

for k=1:numel(I)
    b(:,I(k),J(k))=a(:,k);
end
for k=1:numel(I)
    b(:,J(k),I(k))=a(:,k);
end