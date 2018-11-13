function nt_bsplot(x,w,style,abscissa,zeroflag,rmsflag)
%nt_bsplot(x,sds,style,abscissa,zeroflag,rmsflag) - plot average with bootstrap standard deviation
%
%  x: data to plot (time * trials, or time * 1 * trials)
%  w: weights
%  style: 'zerobased' (default) or 'meanbased'
%  abscissa: use this vector as plot's abscissa (as in 'plot(abscissa,x)' )
%  zeroflag: if 1 draw zero line (default: 1)
%  rmsflag: if 1 use RMS instead of mean (default==0)
%
%  Bootstrap uses N=1000 iterations.
% 
%  Example:
%    nt_bsplot(x)
%  where x is time*trials will plot the average of x over trials, together
%  with +/- 2SDs of the bootstrap resampling.
%
% NoiseTools.
nt_greetings;

if nargin<6 || isempty(rmsflag) ; rmsflag=0; end
if nargin<5 || isempty(zeroflag) ; zeroflag=1; end
if nargin<4; abscissa=[]; end
if nargin<3 || isempty(style); style='zerobased'; end
if nargin<2 w=[]; end
BAND=2; % number of SD in band

x=squeeze(x);
if ~isempty (w); w=squeeze(w); end
if ndims(x)>2; error('X should have at most 2 non-singleton dimensions'); end
[m,n]=size(x);
if n<2; error('bootstrap resampling requires more than 1 column'); end
if isempty(abscissa); abscissa=1:m; end
if numel(abscissa) ~= size(x,1); error('abscissa should be same size as x'); end

N=1000;
if rmsflag
    [a,b]=nt_bsrms(x,N,w);
else
    [a,b]=nt_bsmean(x,N,w);
end
b=b*BAND;


if strcmp(style,'zerobased');
    Y=[b;-flipud(b)]';
elseif strcmp(style,'meanbased');
    Y=[b+a;flipud(-b+a)]';
else
    error('!');
end
abscissa=abscissa(:);
X=[abscissa;flipud(abscissa)];
C=0.7*[1 1 1];
fill(X,Y,C,'edgecolor','none');
hold on;
plot(abscissa,a*0,'k');
plot(abscissa,a, 'b'); 
hold off

% return
% 
% abscissa2=linspace(min(abscissa),max(abscissa),m*2);
% plot(abscissa2,b,'g');
% c=get(gca,'children'); set(c(1),'color',[.7 .7 .7])
% hold on;
% plot(abscissa,a,'b');
% if zeroflag; 
%     plot(abscissa,0*a,'k'); 
%     c=get(gca,'children'); set(c(1),'color',[.5 .5 .5])
% end
% %set(gca,'xlim',[1 m])
% hold off
