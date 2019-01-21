function [iBad,toGood]=nt_find_bad_channels(x,proportion,thresh1,thresh2,thresh3)
%[iBad,toGood]=nt_find_bad_channels(x,proportion,thresh1,thresh2,thresh3) - find bad channels
%
%   iBad: indices of bad channels
%   toGood: matrix to good channels
%
%   x: data
%   proportion: proportion of time above threshold(s) [default: 0.5]
%   thresh1: threshold relative to median absolute value over all data [default: 3]
%   thresh2: absolute threshold
%   thresh3: applies to projection residual
%
% Thresholds apply to absolute value.
% 'thresh3' applies to the residual of the projection of a channel on
% neighboring channels (as calculated by 'sns'), expressed as a proportion
% of median absolute value of that channel. 
%
% NoiseTools


if nargin<2||isempty(proportion); proportion=0.5; end
if nargin<3||isempty(thresh1); thresh1=3; end
if nargin<4; thresh2=[]; end
if nargin<5; thresh3=[]; end

w=ones(size(x));

if ~isempty(thresh1)
    md=median(abs(x(:)));
    w(find(abs(x)>thresh1*md))=0;
end
if ~isempty(thresh2)
    w(find(abs(x)>thresh2))=0;
end
if ~isempty(thresh3)
    NNEIGHBORS=10;
    xx=nt_sns(x,NNEIGHBORS); % remove sensor specific
    md=median(abs(xx(:)));
    xx=xx-x; % residual
    for iChan=1:size(x,2)
        w(find(abs(xx(:,iChan))>thresh3*md),iChan)=0;
    end
end
iBad=find(mean(1-w)>proportion);

toGood=eye(size(x,2));
toGood(:,iBad)=[];

if nargout==0
    % plot, don't return values
    plot(mean(1-w), '.-'); 
    h=line([0 size(x,2)],[proportion proportion]); set(h,'linestyle','--');
    xlabel('channel'); ylabel('proportion bad');
    xlim([0 size(x,2)+1])
    
    clear iBad toGood
end

    
    
    