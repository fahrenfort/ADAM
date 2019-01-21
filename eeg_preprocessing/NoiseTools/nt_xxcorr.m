function [C,idx]=nt_xxcorr(A,B,MAXLAG)
%[C,idx]=nt_xxcorr(A,B,MAXLAG) - true unbiased cross-correlation
%
%  C: unbiased cross-correlation function
%  idx: index of largest extremum.
%  
%  A: first column vector
%  B: second column vector
%  MAXLAG: lags are between -MAXLAG and +MAXLAG.
%

if nargin<1; error('!'); end
if nargin==1; B=A; MAXLAG=[]; end
if nargin==2;
    if numel(B)==1;
        MAXLAG=B; B=A;
    else
        MAXLAG=[];
    end
end
if isempty(MAXLAG); MAXLAG=floor(size(B,1)/4); end
if size(A,1)==1; A=A(:); B=B(:); end
if size(B,1)<=2*MAXLAG; error('!'); end



if size(A,2)==1 && size(B,2)==1; 
    
    % single channels
    B=B(MAXLAG:end-MAXLAG);
    C=xcorr(A,B);
    C=C(size(A,1)-1+(0:2*MAXLAG));
    [~,idx]=max(abs(C));
    
    % power for normalization
    a=cumsum(A.^2);
    a=a(size(B,1)-1:end)-[0;a(1:end-size(B,1)+1)];
    b=sum(B.^2);
    d=sqrt(a*b);
    
    C=C./d; C(find(isnan(C)))=0;

    if nargout==0;
        abscissa=-MAXLAG:MAXLAG;
        plot(abscissa,zeros(size(C)), 'k'); hold on
        plot(abscissa,C);
        plot(abscissa(idx),C(idx),'.k'); hold off
        axis tight; xlabel('lag (samples)');
        lag=idx-MAXLAG-1;
        if C(idx)>1
            if lag>0; 
                title(['X lags Y by ',num2str(lag)]);
            else
                title(['Y lags X by ',num2str(-lag)]);
            end
        else
            if lag>0; 
                title(['X lags -Y by ',num2str(lag)]);
            else
                title(['Y lags -X by ',num2str(-lag)]);
            end
        end
        C=[];
    end
    
else
    
    % multiple channels
    C=zeros(2*MAXLAG+1,size(A,2),size(B,2));
    idx=zeros(size(A,2),size(B,2));
    for k=1:size(A,2)        
        for j=1:size(B,2)
            [a,b]=nt_xxcorr(A(:,k),B(:,j),MAXLAG);
            C(:,k,j)=a;
            idx(k,j)=b;
        end
    end
    
    if nargout==0
       imagescc(idx-MAXLAG-1);
       colorbar
       C=[]; idx=[];
    end
    
end