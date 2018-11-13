function stats=nt_statMatrix(x,plot_params)
%stats=nt_statMatrix(x,plot_params) - calculate statistics arrays for each dim of matrix
%
%  stats: array of statistics arrays
%
%  x: data to stat
%  plot_params: parameters

%% arguments
if nargin<1; error('!'); end
if nargin<2;
    plot_params.bottom=0.1;
    plot_params.height=0.8;
end

%% count true number of dimensions
nDims=ndims(x);
sizeX=size(x);
if nDims==2 && (sizeX(1)==1 || sizeX(2)==1)
    nDims=1;
end
% if we're plotting, concatenate dimensions beyond 4th
if nargout==0 && ndims(x)>4; x=x(:,:,:,:); nDims=4; end


%% statistics on data
stats={}; 
if 0
    % computationally expensive
    for iDim=1:nDims
        sz=size(x);
        x=reshape(x,sz(1),[]);
        stats{iDim}.iDim=iDim;
        stats{iDim}.n=size(x,2);
        stats{iDim}.min=min(x')';
        stats{iDim}.max=max(x')';
        stats{iDim}.mean=nanmean(x')';
        stats{iDim}.var=nanvar(x')';
        x=reshape(x,sz);
        x=shiftdim(x,1);
    end
else
    % cheaper
    sz=size(x);
    for iDim=1:nDims
        stats{iDim}.iDim=iDim;
        stats{iDim}.n=size(x,iDim);
        
        % min
        y=x;
        for k=1:iDim-1;
             y=squeeze(min(y,[],1));
        end
        if iDim<nDims; y=min(y(:,:),[],2); else y=min(y(:,:),[],1); end
        stats{iDim}.min=y(:); 

        % max
        y=x;
        for k=1:iDim-1;
            y=squeeze(max(y,[],1));
        end
        if iDim<nDims; y=max(y(:,:),[],2); else y=max(y(:,:),[],1); end
        stats{iDim}.max=y(:); 

        % mean
        y=x;
        for k=1:iDim-1;
            y=squeeze(nanmean(y,1));
        end
        if iDim<nDims; y=nanmean(y(:,:),2); else y=nanmean(y(:,:),1); end
        stats{iDim}.mean=y(:); 
      
        % var  
        y=x.^2;
        for k=1:iDim-1;
            y=squeeze(nanmean(y,1));
        end
        if iDim<nDims; y=nanmean(y(:,:),2); else y=nanmean(y(:,:),1); end
        stats{iDim}.var=y(:)-stats{iDim}.mean.^2; 
        %nt_whoss
        
    end
end


%% No arguments: put up a window and plot the statistics
if nargout==0;
    if isempty(get(0,'currentfigure'));
        figH=figure;
        set(gcf,'color',[1 1 1]);
    end
    figH=figure(gcf); %clf
    %set(gcf,'name',inputname(1));

    fontsize=12;
    
    % plot in 1,2,3 or 4 panels depending on number of dimensions
    switch nDims
        case 1
            axes('position',[0.05, plot_params.bottom, 0.93, plot_params.height]);  set(gca,'box','on','fontsize',fontsize);
            plot(x, 'k');  title(['n = ', num2str(sizeX(1))]);
            xlabel('samples');
        case 2
            axes('position',[0.05, plot_params.bottom, 0.45, plot_params.height]); set(gca,'box','on','fontsize',fontsize);
            nt_plotstats(stats{1});
            xlabel('samples'); title(['n = ', num2str(sizeX(1))]);
            axes('position',[0.53, plot_params.bottom, 0.45, plot_params.height]); set(gca,'box','on','fontsize',fontsize);
            nt_plotstats(stats{2});
            xlabel('samples'); title(['n = ', num2str(sizeX(2))]); set(gca,'ytick',[]);
        case 3
            axes('position',[0.05, plot_params.bottom, 0.3, plot_params.height]); set(gca,'box','on','fontsize',fontsize);
            nt_plotstats(stats{1});
            xlabel('samples'); title(['n = ', num2str(sizeX(1))]);
            axes('position',[0.37, plot_params.bottom, 0.3, plot_params.height]); set(gca,'box','on','fontsize',fontsize);
            nt_plotstats(stats{2});
            xlabel('samples'); title(['n = ', num2str(sizeX(2))]); set(gca,'ytick',[]);
            axes('position',[0.69, plot_params.bottom, 0.3, plot_params.height]); set(gca,'box','on','fontsize',fontsize);
            nt_plotstats(stats{3});
            xlabel('samples'); title(['n = ', num2str(sizeX(3))]); set(gca,'ytick',[]);
        otherwise % limit to 4 panels (last dims are concatenated)
            axes('position',[0.05, plot_params.bottom, 0.2, plot_params.height]); set(gca,'box','on','fontsize',fontsize);
            nt_plotstats(stats{1});
            xlabel('samples'); title(['n = ', num2str(sizeX(1))]);
            axes('position',[0.27, plot_params.bottom, 0.2, plot_params.height]); set(gca,'box','on','fontsize',fontsize);
            nt_plotstats(stats{2});
            xlabel('samples'); title(['n = ', num2str(sizeX(2))]); set(gca,'ytick',[]);
            axes('position',[0.49, plot_params.bottom, 0.2, plot_params.height]); set(gca,'box','on','fontsize',fontsize);
            nt_plotstats(stats{3});
            xlabel('samples'); title(['n = ', num2str(sizeX(3))]); set(gca,'ytick',[]);
            axes('position',[0.71, plot_params.bottom, 0.2, plot_params.height]); set(gca,'box','on','fontsize',fontsize);
            nt_plotstats(stats{4});
            xlabel('samples'); 
            title_string=['n = ', num2str(sizeX(4))];
            for k=5:numel(sizeX)
                title_string=[title_string,'*',num2str(sizeX(k))];
            end
            title(title_string); 
            set(gca,'ytick',[]);
            
    end
end
            

%%
function nt_plotstats(stats)
%h=nt_plotstats(stats)
%  
%  stats: stats stucture

holdStatus=ishold;
hold on
if isfield(stats,'min')
    nt_plot2(1:size(stats.min,1), [stats.min,stats.max], [1 1 1]*0.9);
end
if isfield(stats,'var')
    sd=sqrt(stats.var);
    nt_plot2(1:size(stats.mean,1), [stats.mean+sd,stats.mean-sd], [1 1 1]* 0.5);
end
nt_plot2(1:size(stats.mean,1),[stats.mean,stats.mean], [1 0 0]);
stats.mean(find(stats.min~=stats.max))=nan;
plot(1:size(stats.mean,1),[stats.mean,stats.mean], '.b');
if holdStatus;
    hold on;
else
    hold off
end

%%
function h=nt_plot2(x,y,c)
%h=nt_plot2(x,y,c) - color region between two plots
%
% 
%  x: abscissa
%  y: ordinate (1 or 2 columns)
%  c: color (see 'fill')
%
%  h: vector of handles to patches
%
%  nt_plot2(x

if nargin<1; error('!'); end

% process parameters
if nargin==1;
    y=x;
    x=(1:size(y,1))';
    c=[0 0 0];
elseif nargin==2;
    c=[0,0,0];
elseif nargin==3;
    ;
else
    error('!');
end

% format data
if size(y,2)==1
    y=[y,zeros(size(y))];
elseif size(y,2)>2
    error('!');
else
    ;
end
x=x(:);

if 0
% make sure that y(:,1)<y(:,2);
    yy=y;
    yy(:,1)=min(y,[],2);
    yy(:,2)=max(y,[],2);
    y=yy;

    % downsample if data array is too large
    TARGET_N=2000;
    if size(x,1)>TARGET_N
        DSR=ceil(size(x,1)/TARGET_N);
        n=floor(size(x,1)/DSR);
        x=x(1:DSR:end);
        x=x(1:n);
        y_extra=y(DSR*n:end,:);
        y=y(1:DSR*n,:);
        a=min(reshape(y(:,1),DSR,size(y,1)/DSR));
        b=max(reshape(y(:,2),DSR,size(y,1)/DSR));
        a(end)=min(a(end),min(y_extra(:,1)));
        b(end)=max(b(end),max(y_extra(:,2)));
        y=[a',b'];
    end
end

% draw plot
yy=flipud(y(:,2));
yy=yy+0.000001*max(abs(y(:))); % workaround for fill bug
h=fill([x;flipud(x)],[y(:,1);yy],c, 'edgecolor', c);
if x(end)>x(1); 
    xlim([x(1)-1,x(end)+1]);
end
set(gca,'box','on');

%% from stats toolbox
function m = nanmean(x,dim)
%% 

% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end

