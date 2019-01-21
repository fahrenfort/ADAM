function [y,stepList]=nt_destep(x,thresh,guard,depth,minstep);
%[y,stepList]=nt_destep(x,thresh,guard,depth,minstep) - remove step glitch from MEG data
%
%  y: step-removed data
%  stepList: indices of steps
%
%  x: data to clean (time * channels)
%  thresh: threshold for variance reduction [default: 0.1]
%  guard: minimum duration of stable interval in samples [default: 1000]
%  depth: recursion depth for nt_split [default:6], determines number of steps
%  minstep: if step size smaller skip step
%
% Searches for large steps, removes them if variance ratio less than 
% "thresh", and both intervals shorter than "guard".
%
% See also nt_detrend, nt_deglitch, nt_split.
%
% NoiseTools.
nt_greetings;

if nargin<1; error('!'); end
if nargin<2||isempty(thresh); thresh=0.1; end
if nargin<3||isempty(guard); guard=100; end
if nargin<4 || isempty(depth); depth=6; end
if nargin<5 ; minstep=[]; end

if isempty(minstep); minstep=(max(x(:))-min(x(:)))*0.0000001; end
y=x;

disp(numel(thresh))

for iChan=1:size(x,2)
    % find step indices
    [stepList,score_vector,score]=nt_split(x(:,iChan),depth,thresh,guard,minstep);
    
    if ~isempty(stepList)
        stepList=[0,stepList,size(x,1)];
        for iSplit=2:numel(stepList)-1
            y1=y(stepList(iSplit-1)+1:stepList(iSplit),iChan); % plateau before
            y2=y(stepList(iSplit)+1:stepList(iSplit+1),iChan); % plateau after
            step=(mean(y2)-mean(y1));
            y(stepList(iSplit)+1:end,iChan)=y(stepList(iSplit)+1:end,iChan)-step;
        end
    end
end
stepList=stepList(2:end-1);


if ~nargout
    % don't return values, just plot
    disp(['steps at: ', num2str(stepList(2:end-1))]);
    figure(1); clf; nt_banner('nt_destep');
    xx=[x(:),y(:)];
    lim=[min(xx(:)),max(xx(:))]; lim=[lim(1)-0.1*(lim(2)-lim(1)), lim(2)+0.1*(lim(2)-lim(1))];
    subplot 211; plot([x,y]); xlabel('samples'); ylim(lim); legend('raw','clean'); legend boxoff
    subplot 212;     plot(y,'r'); xlabel('samples'); legend('clean'); legend boxoff
    clear y, stepList;
end



% test code
% if 0
%     N=8;
%     Wn=0.25; % CTF corner frequency is 1/4 of sampling rate
%     nsamples=100;
%     [B,A] = butter(N,Wn);
%     y=filter(B,A,ones(nsamples,1));
%     BETA0=[1, zeros(1,8), 1, zeros(1,8)];
%     fun = @(beta,x)(filter([beta(1),beta(2),beta(3),beta(4),beta(5),beta(6),beta(7),beta(8),beta(9)],...
%         [beta(10),beta(11),beta(12),beta(13),beta(14),beta(15),beta(16),beta(17),beta(18)],x));
%     x=ones(nsamples,1);
%     BETA = nlinfit(x,y,fun,BETA0);
%     BETA=BETA/BETA(10);
%     BB=BETA(1:9); AA=BETA(10:end);
%     figure(1); clf;
%     plot([y,filter(BB,AA,ones(size(y)))])
% end
if 0
    N=8;
    Wn=0.25; % CTF corner frequency is 1/4 of sampling rate
    [B,A] = butter(N,Wn);
    x=[zeros(1100,1);ones(2000,1)];
    x=filter(B,A,x);
    x=x+0.01*randn(size(x));
    [y,stepList]=nt_destep(x);
    figure(1); clf;
    subplot 211;
    plot([x,y]); legend('raw','destep');legend boxoff
    subplot 212; 
    plot([y,nt_deglitch(y,stepList)]); 
    legend('destep','deglitch');legend boxoff
end
if 0
    N=8;
    Wn=0.25; % CTF corner frequency is 1/4 of sampling rate
    [B,A] = butter(N,Wn);
    x=[zeros(1000,1);ones(1000,1);2*ones(1000,1);3*ones(1000,1);4*ones(1000,1);0*ones(1000,1)];
    x=filter(B,A,x);
    x=x+0.0*randn(size(x));
    [y,stepList]=nt_destep(x);
    figure(1); clf;
    subplot 211;
    plot([x,y]); legend('raw','destep'); legend boxoff
    subplot 212; 
    plot([y,nt_deglitch(y,stepList)]); 
    legend('destep','deglitch');legend boxoff
end
if 0
    x=ft_read_data('/data/meg/litvak/phantom090715_BrainampDBS_20150709_18.ds'); x=permute(x,[2,1,3]); x=x(:,33:306,:);
    x=nt_demean(x(1:2300,100,1));
    x=x*10^12; % -->pT
    [y,stepList]=nt_destep(x);
    figure(1); clf;
    subplot 211;
    plot([x,y]); legend('raw','destep'); legend boxoff
    subplot 212; 
    plot([y,nt_deglitch(y,stepList)]); 
    legend('destep','deglitch');legend boxoff
end

    
