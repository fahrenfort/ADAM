function varargout=nt_spect_plot2(x,varargin)
%nt_spect_plot2 - plot power spectrum
%
%  The power spectral densities of all columns are calculated.  
%  The result is either plotted as an image or returned.
% 
%  See nt_spect_plot, and pwelch for syntax.

if ndims(x)>2;
    [m,n,o]=size(x);
    x=reshape(x,[m,n*o]);
    [pxx,f]=nt_spect_plot2(x,varargin{:});
    pxx=reshape(pxx,[size(pxx,1),n,o]);
    pxx=mean(pxx,3);
else

% ndims(x)==2

    [m,n]=size(x);
    [a,f]=pwelch(x(:,1),varargin{:});
    pxx=zeros(size(a,1),n);
    for k=1:n
        [a,f]=pwelch(x(:,k),varargin{:});
        pxx(:,k)=a;
    end
end

if nargout == 0; % plot
    
    % hack to get nice frequency axis
    [X,f]=nt_spect_plot(x(:,1,1), varargin{:}); % to get scaling factor
    scaling_factor=size(X,1)/max(f);
    nt_spect_plot(x(:,1,1), varargin{:}); % to get x axis labels
    xtick=get(gca,'xtick'); xticklabel=get(gca,'xticklabel');  
    if 0
        pxx=nt_normcol(pxx);
    else
        pxx=bsxfun(@times,pxx,1./max(pxx));
    end
    nt_imagescc(pxx'.^0.25);
    set(gca,'xtick',xtick*scaling_factor,'xticklabel',xticklabel);
    xlabel('Hz'); ylabel('channel');
    
else
    varargout={pxx,f};
end
