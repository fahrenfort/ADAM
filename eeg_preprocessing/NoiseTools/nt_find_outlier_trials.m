function [idx,d]=nt_find_outlier_trials(x,criterion,disp_flag,regress_flag)
%[idx,d,mn,idx_unsorted]=nt_find_outlier_trials(x,criterion,plot,regress_flag) - find outlier trials
%
%  idx: indices of trials to keep
%  d: relative deviations from mean
%  
%  x: data (time * channels * trials)
%  criterion: keep trials less than criterion from mean
%  disp: if true plot trial deviations before and after 
%  regress_flag: if true regress out mean, rather than subtract
%
%  For example criterion=2 rejects trials that deviate from the mean by
%  more than twice the average deviation from the mean.
%

if nargin<2; criterion=inf; end
if nargin<3; disp_flag=1; end
if nargin<4; regress_flag=0; end
if ndims(x)>3; error('x should be 2D or 3D'); end

if ndims(x)==3;
    [m,n,o]=size(x);
    x=reshape(x,m*n,o);
else
    [~,o]=size(x);
end

mn=mean(x,2);
if regress_flag
    mn=nt_tsregress(x,mean(x,2));  % regression
else
    mn=repmat(mn(:),1,o);       % mean
end
d=x-mn; % difference from mean
dd=zeros(1,size(d,2));
for k=1:size(d,2); dd(k)=sum(d(:,k).^2); end
d=dd; clear dd;
d=d/(sum(x(:).^2)/o);

idx=find(d<criterion(1));

if nargout==0;
    % just plot deviations
    plot(d,'.-');
    xlabel('trial'); ylabel('normalized deviation from mean'); 
    clear idx d mn idx_unsorted
else
    if disp_flag
        % plot deviations before & after outlier removal
        figure(100); clf
        nt_banner('outlier trials');
        
        subplot 121; 
        plot(d,'.-'); hold on; 
        plot(setdiff(1:o,idx), d(setdiff(1:o,idx)), '.r');
        xlabel('trial'); ylabel('normalized deviation from mean'); title(['before, ',num2str(numel(d))]);
        drawnow
        
        subplot 122; 
        [~,dd]=nt_find_outlier_trials(x(:,idx),0,[]);
        plot(dd,'.-');
        xlabel('trial'); ylabel('normalized deviation from mean'); title(['after, ',num2str(numel(idx))]);
        drawnow
    end
    
end
     
criterion=criterion(2:end);
if ~isempty(criterion)
    idx=nt_find_outlier_trials(x(:,idx),criterion,disp_flag,regress_flag);
    idx = idx(idx2); % otherwise, idx doesn?t correspond to original matrix anymore
end
