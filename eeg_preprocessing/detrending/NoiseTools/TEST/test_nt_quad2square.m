reset(RandStream.getGlobalStream); % to get reproducible signals

%{
The data are 10-channel data organized into trials.
The stimulus is a single source consisting of a modulated sinusoid with
random phase. It is mixed into the 10 channels via a 1x10 random mixing matrix.
The noise is produced by 9 independent sources mixed into the 10 channels via a
9x10 random mixing matrix.
%}

sr=1000;
nsamples=1000;
ntrials=100;
nchans=10;

% sinusoidal pulse with random phase 
CF1=4;
target=sin(2*pi*0.5*(1:nsamples)'/sr).^2;
target=repmat(target,[1,1,ntrials]) .* ...
    sin( 2*pi*(CF1*repmat((1:nsamples)',[1,1,ntrials])/sr + ...
    repmat(rand(1,1,ntrials),[nsamples,1,1]))); % phase
target=nt_normcol(target);

% noise
NNOISE=9;
noise=randn(nsamples,NNOISE,ntrials);
noise=nt_mmat(noise,randn(NNOISE,nchans));
noise=nt_normcol(noise);

% data
SNR=0.0001;
x=noise+ SNR * nt_mmat(target,randn(1,nchans)) ;
x=nt_demean(x);
x=nt_normcol(x);


%{ 
We append a DC channel (non-zero constant value) and then we form all
cross-products of channels two-by-two.  Appending a DC channel implies
that the set of cross products also includes the original data channels, 
as well as the DC channel.
%}

% append DC channel
x=[x, ones(nsamples,1,ntrials)*mean(abs(x(:)))];

% not sure why this makes things better:
if 1
    THRESH=0;%10^-12;
    [topcs,pwr]=nt_pca0(x,[],[],THRESH); 
    x=nt_mmat(x,topcs);
    x=nt_normcol(x);
end

% form cross products
xx=nt_xprod(x);

% DSS finds the most reproducible quadratic form
keep1=[]; keep2=0;%10^-12;
[toquads,pwr0,pwr1]=nt_dss1(xx,[],keep1,keep2);

% illustrate
figure(1); clf
plot(pwr1./pwr0,'.-')
z=nt_mmat(xx,toquads);
figure(2); clf
subplot 211; nt_bsplot(z(:,1,:));  title('DC')
subplot 212; nt_bsplot(z(:,2,:)); title('most reproducible quadratic form');


%{
We now find the linear component (weighted sum of channels) with square
closest to our optimal quadratic form.
%}

[tosquares,D]=nt_quad2square(toquads(:,2));
z=nt_mmat(x,tosquares);

figure(3); clf
subplot 211; nt_bsplot(z(:,1,:).^2); 
ylabel('power'); title('closest square');
subplot 212; nt_bsplot(z(:,2,:).^2); 
ylabel('power'); title('second closest'); 
figure(4); clf
plot(abs(D(2:end)))

%{
It may necessary to sort the components to ensure that the DC component
does not come first.
%}

