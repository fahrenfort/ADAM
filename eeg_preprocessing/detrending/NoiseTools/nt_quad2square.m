function [tosquares,D]=nt_quad2square(toquad,order)
%[tosquare,D]=nt_quad2square(toquad,order) - quadratic to squared linear component
% 
%   tosquares: matrix to linear components with closest square
%   D: scores
%
%   toquad: matrix that defines quadratic component from cross products
%   order: order of cross products, 'colwise', 'diagwise' (default)

if nargin<2; order='diagwise'; end
if nargin<1; error('!'); end

if size(toquad,2)~=1; error('toquad should be column vector'); end

nquads=size(toquad,1);
nchans=floor(sqrt(nquads*2));
if nquads ~= nchans*(nchans+1)/2; 
    [nchans nquads]
    error('unexpected size for toquad');
end

ii=1;
A=zeros(nchans);
switch order
    case 'colwise'
        for iRow=1:nchans;
            for iCol=1:iRow
                if iRow==iCol;
                    A(iRow,iCol)=toquad(ii);
                else
                    A(iRow,iCol)=toquad(ii)/2;
                    A(iCol,iRow)=toquad(ii)/2;
                end
                ii=ii+1;
            end
        end
    case 'diagwise'
        for iDiag=0:nchans-1
            for kk=1:(nchans-iDiag)
                iRow=kk+iDiag;
                iCol=kk;
                if iRow==iCol;
                    A(iRow,iCol)=toquad(ii);
                else
                    A(iRow,iCol)=toquad(ii)/2;
                    A(iCol,iRow)=toquad(ii)/2;
                end
                ii=ii+1;
            end
        end
    otherwise
        error('unexpected order');
end

% eigenvectors & values
[V,D]=eig(A);
D=diag(D);
[~,idx]=sort(abs(D),'descend'); 
D=D.^2/sum(D.^2);
V=V(:,idx);
D=D(idx);
tosquares=V;
