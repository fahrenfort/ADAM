function nt_bsplot2(x,percentile,style,abscissa,zeroflag,rmsflag)
%nt_bsplot(x,percentile,style,abscissa,zeroflag,rmsflag) - plot average with confidence interval
%
%  x: data to plot (time * trials, or time * 1 * trials)
%  percentile: width of band to plot in standard deviations (default: 2)
%  style: 'zerobased' (default) or 'meanbased'
%  abscissa: use this vector as plot's abscissa (as in 'plot(abscissa,x)' )
%  zeroflag: if 1 draw zero line (default: 1)
%  rmsflag: if 1 use RMS instead of mean (default==0)
% 
%  Example:
%    nt_bsplot(x)
%  where x is time*trials will plot the average of x over trials, together
%  with +/- 2SDs of the bootstrap resampling.
%
% NoiseTools.

if nargin<6 || isempty(rmsflag) ; rmsflag=0; end
if nargin<5 || isempty(zeroflag) ; zeroflag=1; end
if nargin<4; abscissa=[]; end
if nargin<3 || isempty(style); style='zerobased'; end
if nargin<2 || isempty(percentile); percentile=[5 95]; end

if numel(percentile)==1
    percentile=[percentile,100-percentile];
end

x=squeeze(x);
if ndims(x)>2; error('X should have at most 2 non-singleton dimensions'); end
[m,n]=size(x);
if n<2; error('bootstrap resampling requires more than 1 column'); end
if isempty(abscissa); abscissa=1:m; end
if numel(abscissa) ~= size(x,1); error('abscissa should be same size as x'); end

N=1000;
if rmsflag
    [a,~,c]=nt_bsrms(x,N);
else
    [a,~,c]=nt_bsmean(x,N);
end
percentile
c=[prctile(c',percentile(1))',prctile(c',percentile(2))'];
if strcmp(style,'zerobased')
    c=c-repmat(a,1,2);
elseif strcmp(style,'meanbased');
    ;
else
    error('!');
end
Y=[c(:,1);flipud(c(:,2))];
    
abscissa=abscissa(:);
X=[abscissa;flipud(abscissa)];
C=0.7*[1 1 1];
fill(X,Y,C,'edgecolor','none');
hold on;
plot(a*0,'k');
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
