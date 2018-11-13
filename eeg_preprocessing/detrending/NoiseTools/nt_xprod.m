function [y,ind]=nt_xprod(x,flag,dsratio,normrow_flag)
%[y,ind]=nt_xprod(x,flag,dsratio,normrow_flag) - form all crossproducts
%
%  y: crossproducts 
%  ind: linear index of cross-products 
% 
%  x: data (time*channels*trials)
%  flag: 'lower','nodiag','full' [default: 'lower']
%  dsratio: ratio by which to downsample cross-product [default: 1]
%  normrow_flag: if true, divide each slice by trace (default: false)
%
% If flag is 'lower' (default), return lower diagonal terms, including
% diagonal.  The order of cross-products is diagonal first (squares first). 
% 
% If 'nodiag' return lower diagonal terms without diagonal.
%
% If 'full', return full array of cross-products as a 3D matrix.
%
%

if nargin<4 || isempty(normrow_flag); normrow_flag=0; end
if nargin<3 || isempty(dsratio); dsratio=1; end
if nargin<2 || isempty(flag); flag='lower'; end

if ~strcmp(flag,{'lower','nodiag','full'}); error('unexpected flag'); end

if ndims(x)==3
    if rem(size(x,1),dsratio); error('size(x,1) must be multiple of dsratio'); end
    y=nt_fold(nt_xprod(nt_unfold(x),flag,dsratio), size(x,1)/dsratio);
else
    [nsamples,nchans]=size(x);
    nsamples=floor(nsamples/dsratio);

    y=zeros(nsamples,nchans*(nchans+1)/2);
    ind=zeros(nchans*(nchans+1)/2,1);
    start=0;
    
    iProd=1;
    for iDiag=start:nchans-1
        for kk=1:(nchans-iDiag)
            xx=x(:,kk+iDiag).*x(:,kk);
            y(:,iProd)=nt_dsample(xx,dsratio);
            ind(iProd)=sub2ind([nchans,nchans],kk+iDiag,kk);
            iProd=iProd+1;
        end
    end
    if normrow_flag
        for iRow=1:size(y,1)
            y(iRow,:)=y(iRow,:)/(eps+sum(y(iRow,1:nchans)));
        end
    end
    if strcmp(flag, 'nodiag')
        y=y(:,nchans+1:end);
    end

    % this could be optimized to save memory:
    if strcmp(flag,'full')
        y0=y;
        y=zeros(nsamples,nchans,nchans);

        [I,J]=ind2sub(nchans,ind);

        for k=1:numel(I)
            y(:,I(k),J(k))=y0(:,k);
            y(:,J(k),I(k))=y0(:,k);
        end

        ind=[];
    end
    
%     switch order
%         case 'colwise'
%             for iRow=1:nchans
%                 for iCol=1:iRow
%                     xx=x(:,iCol).*x(:,iRow);
%                     y(:,iProd)=nt_dsample(xx,dsratio);
%                     ind(iProd)=sub2ind([size(x,2),size(x,2)],iRow,iCol);
%                     iProd=iProd+1;
%                 end
%             end
%         case 'diagwise'
%             for iDiag=start:nchans-1
%                 for kk=1:(nchans-iDiag)
%                     xx=x(:,kk+iDiag).*x(:,kk);
%                     y(:,iProd)=nt_dsample(xx,dsratio);
%                     ind(iProd)=sub2ind([size(x,2),size(x,2)],kk+iDiag,kk);
%                     iProd=iProd+1;
%                 end
%             end
%         otherwise
%             error('unexpected order flag');
%     end
end


