function [d,err_rate,perm_rate,area_roc,roc]=nt_dprime(x,y,jd_flag)
%[d,err_rate,perm,area_roc,roc]=nt_dprime(x,y,jd_flag) - calculate d' (discriminability) of two distributions
%
%  d: discriminablity index
%  err_rate: error rate for linear discrimination
%  perm_rate: rate for permuted data
%  area_roc: area under ROC curve
%  roc: ROC curve
%
%  x, y: data (column vectors or matrices)
%  jd_flag: apply JD first
% 
% NoiseTools

NSTEPS=1000; % number of steps to find min of error rate
P=0.05; % theshold value for permutation test
NPERMUTE=1000; % number of trials for permutation test

if nargin<3; jd_flag=[]; end
if nargin<2; error('!'); end
if iscell(x)
    xx=[]; yy=[];
    for iCell=1:numel(x);
        xx=[xx; x{iCell}];
        yy=[yy; y{iCell}];
    end
    x=xx; y=yy;
end
if size(x,2) ~= size(y,2); error('!'); end

if jd_flag; 
    c0=nt_cov(nt_demean(x))+nt_cov(nt_demean(y));
    c1=nt_cov(mean(x)-mean(y));
    todss=nt_dss0(c0,c1);
    x=nt_mmat(x,todss);
    y=nt_mmat(y,todss);
end

d=abs(mean(x)-mean(y)) ./ sqrt((var(x)+var(y))/2);

% make sure that y>x
for iChan=1:size(x,2)
    if mean(x(:,iChan))>mean(y(:,iChan));
        x(:,iChan)=-x(:,iChan);
        y(:,iChan)=-y(:,iChan);
    end
end

if nargout>1; % error rate
    err_rate=[];
    for iChan=1:size(x,2)
        min_error=1;
        for thresh=linspace(mean(x(:,iChan)),mean(y(:,iChan)),NSTEPS);
            x2y=sum(x(:,iChan)>thresh);
            y2x=sum(y(:,iChan)<thresh);
            nErr=x2y+y2x;
            min_error=min(min_error,nErr/(size(x,1)+size(y,1)));
        end
        err_rate(iChan)=min_error;
    end
end

if nargout>2; %permutation test
    perm_rate=[];
    for iChan=1:size(x,2)
        for iPermute=1:NPERMUTE
            if rem(iPermute,100)==1; disp(iPermute); end
            % scramble between x and y:
            z=[x(:,iChan);y(:,iChan)]; 
            z=z(randperm(size(z,1)));
            xx=z(1:size(x,1));    
            yy=z(size(x,1)+1:end);
            min_error=1;
            % scan criterion for minimum error
            for thresh=linspace(mean(x(:,iChan)),mean(y(:,iChan)),NSTEPS);
                x2y=sum(xx>thresh);
                y2x=sum(yy<thresh);
                nErr=x2y+y2x;
                min_error=min(min_error,nErr/(size(x,1)+size(y,1)));
            end
            min_errors(iPermute)=min_error;
        end
        % find 5th percentile of distribution of error rates
        min_errors=sort(min_errors);
        min_errors=min_errors(1:round(NPERMUTE*P));
        perm_rate(iChan)=min_errors(end);
    end
end

if nargout>3; error('ROC not implemented yet'); end


% test code
if 0
    x=randn(10000,1);
    y=1+randn(10000,1);
    figure(1); clf
    t=-3:0.1:4;
    plot(t,hist(x,t));
    hold on;
    plot(t,hist(y,t), 'r');
    [d,e,p]=nt_dprime(x,y);
    disp(['d'': ', num2str(d)]);
    disp(['e'': ', num2str(e)]);
    disp(['p'': ', num2str(p)]);
end
