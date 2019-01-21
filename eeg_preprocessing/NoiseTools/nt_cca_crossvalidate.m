function [AA,BB,RR,SD]=nt_cca_crossvalidate(xx,yy,shifts,doSurrogate)
%[AA,BB,RR,SD]=nt_cca_crossvalidate(xx,yy,shifts,doSurrogate) - CCA with cross-validation
%
%  AA, BB: cell arrays of transform matrices
%  RR: r scores (2D)
%  SD: standard deviation of correlation over non-matching pairs (2D)
%
%  xx,yy: cell arrays of column matrices
%  shifts: array of shifts to apply to y relative to x (can be negative)
%  doSurrogate: if true estimate sd of correlation over non-matching pairs

if nargin<4; doSurrogate=[]; end
if nargin<3; shifts=[0]; end
if nargin<2; error('!'); end
if ~iscell(xx) || ~iscell(yy); error('!'); end
if length(xx) ~= length (yy); error('!'); end
if size(xx{1},1) ~= size(yy{1},1); error('!'); end

if nargout==0 || nargout==4; doSurrogate=1; end

%%
% calculate covariance matrices
nTrials=length(xx);
n=size(xx{1},2)+size(yy{1},2);
C=zeros(n,n,length(shifts),nTrials);
disp('Calculate all covariances...');
nt_whoss;
for iTrial=1:nTrials
    C(:,:,:,iTrial)=nt_cov_lags(xx{iTrial}, yy{iTrial},shifts);
end

%%
% calculate leave-one-out CCAs
disp('Calculate CCAs...');
for iOut=1:nTrials
    CC=sum(C(:,:,:,setdiff(1:nTrials,iOut)),4); % covariance of all trials except iOut
    [A,B,R]=nt_cca([],[],[],CC,size(xx{1},2));  % corresponding CCA
    AA{iOut}=A;
    BB{iOut}=B;
end
clear C CC

%%
% calculate leave-one-out correlation coefficients
disp('Calculate cross-correlations...');
for iOut=1:nTrials
    iNext=mod(iOut,nTrials)+1; % correlate with next in list
    A=AA{iOut};
    B=BB{iOut};
    for iShift=1:length(shifts)
        [x,y]=nt_relshift(xx{iOut},yy{iOut},shifts(iShift));
        a=A(:,:,iShift);
        b=B(:,:,iShift);
        r(:,iShift)=diag( nt_normcol(x*a)' * nt_normcol(y*b )) / size(x,1); 
    end
    RR(:,:,iOut)=r;
    if doSurrogate
        for iShift=1:length(shifts)
            [x,y]=nt_relshift(xx{iOut},yy{iNext},shifts(iShift));
            a=A(:,:,iShift);
            b=B(:,:,iShift);
            mn=min(size(x,1),size(y,1));
            s(:,iShift)=diag( nt_normcol(x(1:mn,:)*a)' * nt_normcol(y(1:mn,:)*b )) / mn; 
        end
        ss(:,:,iOut)=s;
    end
end
if doSurrogate
    VAR=(sum(ss.^2,3)-sum(ss,3).^2/nTrials) / (nTrials-1);
    SD(:,:)=sqrt(VAR);
end
disp('done');

%%
% If no output arguments, plot something informative

if nargout==0
    figure(1); clf;
    if length(shifts)>1; 
        plot(mean(RR,3)'); title('correlation for each CC'); xlabel('shift'); ylabel('correlation');
        hold on; 
        plot(SD', ':r');
        legend('correlation','standard error'); legend boxoff
    else
        plot(squeeze(mean(RR,3))); title ('correlation for each CC'); xlabel('CC'); ylabel('correlation');
        plot(SD', ':r');
    end
    figure(2); clf;
    size(RR)
    for k=1:min(4,size(RR,1))
        subplot(2,2,k);
        [~,idx]=max(mean(RR(k,:,:),3));
        [x,y]=nt_relshift(xx{1},yy{1},shifts(idx));
        plot([x*A(:,k,idx), y*B(:,k,idx)]);
        disp(corr(nt_normcol([x*A(:,k,idx), y*B(:,k,idx)])));
        title(['CC ',num2str(k)]); xlabel('sample'); 
    end
end

