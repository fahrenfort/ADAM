function [bc,K]= rst_RASH(varargin)

% Computes the density estimate of data using ASH or RASH
% - Average Shifted Histogram is taken from Martinez and Martinez - Comp
% statistics handbook with Matlab.
% - Random Average Shifted Histogram algorithm is coded based on Bourel et
% al. Computational Statistics and Data Analysis 79 (2014) 149–164
%
% FORMAT [bc,K]=rst_RASH(data)
%        [bc,K]=rst_RASH(data,m,type)
%
% INPUT  data is a vector
%        m is how many histograms to compute (Default = 100);
%        type is 'ASH' or 'RASH'
%
% OUTPUT if none specified it creates a figure
%        bc is the bin count
%        K is the estimated histogram / kernel
%
% Cyril Pernet 2 August 2016
% The University of Edinburgh
% -----------------------------------
% Copyright Robust Statistical Toolbox

data = varargin{1};
if sum(isnan(data)) ~=0
    warndlg('NaN data included')
end

if size(data,1) == 1
    data = data';
end
n = length(data);
m = 100; % this is the parameter setting the number of hist in ASH and RASH
type = 'RASH'; % method to use

if nargin == 2
    m = varargin{2};
elseif nargin == 3
    m = varargin{2};
    type = varargin{3};
end
clear varargin

% % standard probability histogram
% % -----------------------------
% % 1 determine how many bins
% h = 2.15*sqrt(var(data))*n^(-1/5);
% t0 = min(data)-1;
% tm = max(data)+1;
% bins = t0:h:tm;
% % 2 compute frequency
% nu=hist(data,bins);
% nu(end) = [];
% % 3 make it as prob
% nu = nu./(n*h);
% % 4 plot
% bc = (t0+h/2):h:(tm-h/2);
% figure;
% bar(bc,nu,1,'FaceColor',[0.5 0.5 1]);
% title('Prob histogram')

% a typooical problem is with data that have a small range like
% probabilities, the meash is not of the right size or contains
% to few bins - the silution here is simply to scale the data
if range(data) <= 1
    data = data.*10; 
    scale = 1;
else
    scale = 0;
end

switch type
    case('ASH')
        
        % Average Shifted Histogram
        % ----------------------------------
        h = 2.15*sqrt(var(data))*n^(-1/5);
        delta = h/m;
        % 1 make a mesh with size delta
        offset = min(diff(data))/2; 
        if abs(offset) > 1
            offset = 0.5;
        end
        t0 = min(data) - offset;
        tf = max(data) + offset;
        nbin = ceil((tf-t0)/delta);
        binedge = t0:delta:(t0+delta*nbin);
        out = find(binedge>tf);
        if out == 1
            binedge(out) = tf;
        else
            binedge(out(1)) = tf;
            binedge(out(2:end)) = [];
        end
        % 2 bin count
        nu = histc(data,binedge);
        nu = [zeros(1,m-1) nu' zeros(1,m-1)];
        % 3 Get the weight vector.
        kern = inline('(15/16)*(1-x.^2).^2');
        ind = (1-m):(m-1);
        den = sum(kern(ind/m));% Get the denominator.
        wm = m*(kern(ind/m))/den;% Create the weight vector.
        % 4 Get the bin heights over smaller bins.
        ASH=zeros(1,nbin);
        for k=1:nbin
            ind=k:(2*m+k-2);
            ASH(k)=sum(wm.*nu(ind));
        end
        K = ASH/(n*h);
        bc = t0+((1:nbin)-0.5)*delta;
        if scale == 1
            bc = bc./10;
        end
        
        if nargout == 0
            figure; bar(bc,K,1,'FaceColor',[0.5 0.5 1]);
            title(['ASH m = ' num2str(m)]); axis tight;
            grid on; box on;
        end
        
    case('RASH')
        % Random Average Shifted Histogram
        % ----------------------------------
        
        h = 2.15*sqrt(var(data))*n^(-1/5);
        delta = h/m;
        % 1 make a mesh with size delta
        offset = min(diff(data))/2; 
        if abs(offset) > 1
            offset = 0.5;
        end
        t0 = min(data) - offset;
        tf = max(data) + offset;
        nbin = ceil((tf-t0)/delta);
        binedge = t0:delta:(t0+delta*nbin);
        out = find(binedge>tf);
        if ~isempty(out)
            if length(out)== 1
                binedge(out) = tf;
            else
                binedge(out(1)) = tf;
                binedge(out(2:end)) = [];
            end
        end
        % 2 Get the weight vector.
        kern = inline('(15/16)*(1-x.^2).^2');
        ind = (1-m):(m-1);
        den = sum(kern(ind/m));% Get the denominator.
        wm = m*(kern(ind/m))/den;% Create the weight vector.
        % 3 compute bin with shifted edges
        RH=zeros(1,nbin);
        RSH=zeros(m,nbin);
        for e=1:m
            v = binedge + (delta*randn(1,1)); % e is taken from N(0,h);
            v(v<t0) = t0; % lower bound
            v(v>tf) = tf; % upper bound
            nu = histc(data,v);
            nu = [zeros(1,m-1) nu' zeros(1,m-1)];
            for k=1:nbin
                ind=k:(2*m+k-2);
                RH(k)=sum(wm.*nu(ind));
            end
            RSH(e,:) = RH/(n*h);
        end
        K = mean(RSH,1);
        bc = t0+((1:nbin)-0.5)*delta;
        offset = min(bc)-min(data);
        bc = bc - offset;
        if scale == 1
            bc = bc./10;
        end
        
        if nargout == 0
            figure; bar(bc,K,1,'FaceColor',[0.5 0.5 1]);
            title(['RASH m = ' num2str(m)]);
            grid on; box on;
        end

end
