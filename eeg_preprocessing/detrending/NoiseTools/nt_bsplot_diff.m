function nt_bsplot_diff(x,y,band,style,abscissa,zeroflag)
%nt_bsplot_diff(x,y,sds,style,abscissa,zeroflag,rmsflag) - plot average difference with bootstrap standard deviation
%
%  x, y: data (time * trials, or time * 1 * trials)
%  band: width of band to plot in standard deviations (default: 2)
%  style: 'zerobased' (default) or 'meanbased'
%  abscissa: use this vector as plot's abscissa (as in 'plot(abscissa,x)' )
%  zeroflag: if 1 draw zero line (default: 1)
%
%  Bootstrap uses N=1000 iterations.
% 
%  Example:
%    nt_bsplot_diff(x,y)
%  where x and y are time*trials will plot the average of x - y over trials, together
%  with +/- 2SDs of the bootstrap resampling.
%
% NoiseTools.

if nargin<6 || isempty(zeroflag) ; zeroflag=1; end
if nargin<5; abscissa=[]; end
if nargin<4 || isempty(style); style='zerobased'; end
if nargin<3 || isempty(band); band=2; end
if nargin<2; error('!'); end

x=squeeze(x);
if ndims(x)>2; error('X should have at most 2 non-singleton dimensions'); end
[m,n]=size(x);
if n<2; error('bootstrap resampling requires more than 1 column'); end
y=squeeze(y);
if ndims(y)>2; error('Y should have at most 2 non-singleton dimensions'); end
[m2,n2]=size(y);
if n2<2; error('bootstrap resampling requires more than 1 column'); end

if m ~= m2; error('X and Y should have same number of rows'); end
    
if isempty(abscissa); abscissa=1:m; end
if numel(abscissa) ~= size(x,1); error('abscissa should be same size as first dim of x and y'); end

N=1000;
[a,b]=nt_bsmean_diff(x,y,N);


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
plot(abscissa,a); 
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
