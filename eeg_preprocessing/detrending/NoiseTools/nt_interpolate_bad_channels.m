function [toGood,fromGood]=interpolate_bad_channels(x,iBad,coordinates,n)
%y=interpolate_bad_channels(x,iBad,coordinates,n) - interpolate bad channels from good
%
%  y: interpolated data
% 
%  x: data to interpolate
%  iBad: indices of bad channels
%  coordinates: coordinate map (see nt_proximity)
%  n: number of neighboring channels to use [default: 3]
%
% NoiseTools;

nt_greetings;

if nargin<3; 
    error('!'); 
end
if nargin<4; 
    n=3;
end

nchans=size(x,2);
toGood=eye(nchans);
toGood(:,iBad)=[];

[closest,d]=nt_proximity(coordinates);
if size(closest,1)~=nchans; error('!'); end

fromGood=eye(nchans);
for iChan=iBad
    iOthers=closest(iChan,:);
    iOthers=setdiff(iOthers, iBad, 'stable'); % don't include bad channels
    if numel(iOthers)<n; error('!'); end
    iOthers=iOthers(1:n);
    w=1./(d(iChan,iOthers) + eps);
    w=w/sum(w);
    fromGood(iOthers,iChan)=w;
end
fromGood(iBad,:)=[];
    
topo=ones(nchans,1);
topo(iBad)=0;
if nargout==0
    figure(100); clf
    subplot 121; nt_imagescc(fromGood);
    subplot 122; nt_topoplot(coordinates,topo);
end

