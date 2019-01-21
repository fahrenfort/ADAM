function varargout=nt_spect_plot(x,varargin)
%nt_spect_plot - plot power spectrum
%
%  The power spectral densities of all columns (pages, etc.) are calculated
%  and added.  The result is divided by the number of columns and either
%  plotted or returned.
% 
% See pwelch for syntax.

if numel(x)==0; error('!'); end

N=numel(x);
ncols=N/size(x,1);
x=reshape(x,size(x,1),ncols);

[pxx,f]=pwelch(x(:,1),varargin{:});

for k=1:ncols
    [a,b]=pwelch(x(:,k),varargin{:});
    pxx=pxx+a;
end
pxx=pxx/ncols;


if nargout == 0;
    plot(f,abs(pxx).^0.5);
    set(gca,'yscale','log');
    xlim([f(1) f(end)]);
    xlabel('Hz'); ylabel('Hz^{-0.5}');
    varargout={};
else
    varargout={pxx,f};
end
