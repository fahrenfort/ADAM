function [squares,quads,D]=nt_qca(x,npcs,nsmooth,nquads)
%[squares,quads,D]=nt_qca(x,npcs,nsmooth,nquads) - maximize induced power using quadratic component analysis
%
%  squares: linear components closest to most repeatable quadratic component 
%  quads: most reproducible quadratic component(s)
%  D: eigenvalues
%
%  x: data (time*channel*trial)
%  npcs: maximum number of data PCs to use (if [] use all)
%  nsmooth: square smoothing window to apply to xproducts [default: 1]
%  nquads: number of quadratic components to return [default: 1]
%
%  Usually we are interested in the first component of 'squares'
%  (component with square closest to best quadratic component). 
% 
% See nt_qca0, nt_quad2square.
% NoiseTools.


if nargin<4||isempty(nquads); nquads=1; end
if nargin<3||isempty(nsmooth); nsmooth=1; end
if nargin<2; error('!'); end
[nsamples,nchans,ntrials]=size(x);

x=[x,ones(nsamples,1,ntrials)*max(abs(x(:)))]; % append a DC component to absorb DC

if nargout==1;
    tosquares=nt_qca0(x,npcs,nsmooth,nquads);
else
    [tosquares,quads,D]=nt_qca0(x,npcs,nsmooth,nquads);
    quads=quads(:,2:end,:); % discard first (DC)
end

squares=nt_mmat(x,tosquares);
squares=nt_demean2(squares);

if 0
    r=nt_repeatability(squares.^2);
    [dummy,idx]=sort(r,'descend');
    squares=squares(:,idx,:);
end

if nargout==0;
    disp('no output arguments: plot');
    figure(100); 
    subplot 311; 
    plot(abs(D), '.-'); xlabel('component'); ylabel('score');
    subplot 312; 
    nt_bsplot(quads(:,1,:));
    title('best quadratic');
    subplot 313; 
    nt_bsplot(squares(:,1,:).^2);
    xlabel('samples'); title('closest square');
    clear squares quads D
end





    