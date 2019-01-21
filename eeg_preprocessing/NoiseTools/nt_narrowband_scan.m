function A=nt_narrowband_scan(x,freqs,sr,Q,plotflag,yulewalk_order)
%A=nt_narrowband_scan(x,freqs,sr,Q,plotflag) - scan for narrowband components using DSS
%
%  A: cell array of DSS matrices
%
%  x: data (time*channels or time*channels*trials)
%  freqs: Hz, array of bias frequencies
%  sr: Hz, sampling rate
%  Q: quality factors of scanning filter (default: [8 4])
%  plotflag: if true plot (default: 1)
%  yulewalk_order: if present apply inverse filtering with this order
%
% If no output arguments, plots spectra of first DSS components for each
% bias
%
nt_greetings;

if nargin<6; yulewalk_order=[]; end
if nargin<5; plotflag=[]; end
if nargin<4||isempty(Q); Q=[8 4]; end
if nargin<3; error('!'); end

freqs=freqs(:);


current_figure=get(0,'CurrentFigure');
A={};
for iFreqs=1:numel(freqs)
    freq=freqs(iFreqs);
    [b,a]=nt_filter_peak(freq/(sr/2),Q(1));
    [c0,c1]=nt_bias_filter(x,b,a);
    if numel(Q)==2
        [b,a]=nt_filter_peak(freq/(sr/2),Q(2));
        [~,c0]=nt_bias_filter(x,b,a);
    end
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    
    if 1; figure(100); clf; plot(pwr1./pwr0,'.-'); title([num2str(freq), 'Hz bias']); ylabel('score');  drawnow;end
    
    A{iFreqs}=todss;
end


if ~nargout || ~isempty(plotflag)
    if isempty(current_figure);
        figure; 
    else
        figure(current_figure);
    end
    AA=zeros(size(x,2),numel(freqs));
    for iFreqs=1:numel(freqs)
        AA(:,iFreqs)=A{iFreqs}(:,1);
    end
    x=nt_mmat(x,AA);
    nfft=2.^nextpow2(size(x,1)+1)/2;
    MAX_NFFT=1024;  
    nfft=min(nfft,MAX_NFFT);
    nt_spect_plot2(nt_normcol(x),nfft,[],[],sr);
    K=round(numel(freqs)/6);
    set(gca,'ytick',1:K:numel(freqs), 'yticklabel',num2str(freqs(1:K:end), '%.3g')); ylabel('Hz');
    set(gca,'xgrid','on','xminortick','on');
    drawnow;
end

if ~nargout;     clear A; end
