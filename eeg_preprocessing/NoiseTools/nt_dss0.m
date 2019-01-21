function [todss,pwr0,pwr1]=nt_dss0(c0,c1,keep1,keep2)
%[todss,pwr1,pwr2]=nt_dss0(c0,c1,keep1,keep2) - dss from covariance
%
% todss: matrix to convert data to normalized DSS components
% pwr0: power per component (baseline)
% pwr1: power per component (biased)
%
% c0: baseline covariance
% c1: biased covariance
% keep1: number of PCs to retain (default: all)
% keep2: ignore PCs smaller than keep2 (default: 10.^-9)
%
% NoiseTools
nt_greetings;

if nargin<4||isempty(keep2); keep2=10.^-9; end
if nargin<3; keep1=[]; end
if nargin<2; error('needs at least two arguments'); end

if size(c0)~=size(c1); error('C0 and C1 should have same size'); end
if size(c0,1)~=size(c0,2); error('C0 should be square'); end

if any(find(isnan(c0)))
    error('NaN in c0');
end
if any(find(isnan(c1)))
    error('NaN in c1');
end
if any(find(isinf(c0)))
    error('INF in c0');
end
if any(find(isinf(c1)))
    error('INF in c1');
end
% PCA and whitening matrix from the unbiased covariance
[topcs1,evs1]=nt_pcarot(c0,keep1,keep2);
evs1=abs(evs1);

% truncate PCA series if needed
if ~isempty(keep1); topcs1=topcs1(:,1:keep1); evs1=evs1(1:keep1); end
if ~isempty(keep2); idx=find(evs1/max(evs1)>keep2); topcs1=topcs1(:,idx); evs1=evs1(idx);  end

% apply PCA and whitening to the biased covariance
N=diag(sqrt(1./(evs1)));     
c2=N'*topcs1'*c1*topcs1*N;

% matrix to convert PCA-whitened data to DSS
[topcs2,evs2]=nt_pcarot(c2,keep1,keep2);

% DSS matrix (raw data to normalized DSS)
todss=topcs1*N*topcs2;
N2=diag(todss'*c0*todss);
todss=todss*diag(1./sqrt(N2)); % adjust so that components are normalized


% power per DSS component
pwr0=sqrt(sum((c0'*todss).^2)); % unbiased
pwr1=sqrt(sum((c1'*todss).^2)); % biased



