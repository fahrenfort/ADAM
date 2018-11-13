function y=nt_deboing(x,events)
%y=nt_deboing(x,events) - fit, remove ringing associated with events
%
%  y: clean data
%
%  x: data to clean
%  events: list of event samples
%
% NoiseTools.
nt_greetings;

if nargin<2; error('!'); end

if size(x,2)>1; 
    error('x should be column vector'); 
end

ORDER=10; % order of polynomial trend
NSAMPLES=100; % number of samples over which to estimate impulse response
EXTRA=50; % samples before stimulus to anchor trend
NNUM=8;NDEN=8; % number of filter coeffs
THRESH=3; % threshold for robust detrending

% remove events too close to beginning or end
events(find(events<EXTRA))=[];
events(find(events>size(x,1)-NSAMPLES))=[];

y=x;
for iEvent=1:numel(events)
    event=events(iEvent);
    
    % select portion to fit filter response, remove polynomial trend
    event_response=x(event-EXTRA:event+NSAMPLES);
    event_response=nt_detrend(event_response,ORDER,[],[],THRESH); 
    event_response=event_response(EXTRA+1:end);
       
    event_response=event_response(2:end); % not sure why...
    
    % estimate filter parameters
    event_response=[event_response;zeros(size(event_response))]; % helps ensure stable filter  
    [B,A]=stmcb(event_response,NNUM,NDEN); 
    
    % estimate filter response to event
    model=filter(B,A,[1;zeros(NSAMPLES,1)]);
    y(event+(1:size(model,1)))=x(event+(1:size(model,1))) - model;
    %figure(1); clf; plot(nt_demean([event_response(1:size(model,1)),model,x(event+(1:size(model,1)))])); pause
end

if nargout==0
    % plot, don't return
    figure(1); clf;
    plot([x,y]);
    
    clear y
end

        
