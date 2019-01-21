function [y,mn]=nt_rereference(x,w,factor)
%[y,mn]=nt_rereference(x,thresh,w,factor) - rereference by subtracting robust mean
%
%   y: rereferenced data
%   mn: subtracted reference
%
%   x: data to rereference (time*channels or time*channels*trials
%   w: mask
%   factor: apply to mean before subtracting [default: 1]
%   
% See also: nt_detrend, nt_outliers
%
% NoiseTools

if nargin<1; error('!'); end
if nargin<2||isempty(w); w=ones(size(x)); end
if nargin<3||isempty(factor); factor=1; end

if size(w,1)~=size(x,1); error('!'); end
if size(w,2)==1; w=repmat(w,1,size(x,2)); end
if size(w,2) ~= size(x,2); error('!'); end

%{
Reference is weighted mean.
%}

mn=sum(x.*w,2)./sum(w,2);
mn(isnan(mn))=0;
y=bsxfun(@minus,x,mn);

if ~nargout
    % don't return values, just plot
    figure(1); clf;
    subplot 411; plot(x); title('raw');
    subplot 412; plot(y); title('rereferenced');
    subplot 413; plot(mn); title('reference');
    subplot 414; plot(mean(w)); title('mean weight'); ylim([-.1 1.1']); xlabel('channel');
    figure(2); 
    histogram(mn(:)/std(mn(:)), 100);
    clear y mn
end