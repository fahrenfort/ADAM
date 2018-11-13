clear
N=10000;

% single channel, step
x=[randn(N,1)+1; randn(N,1)-1]; 

figure(1); clf;
subplot 411
plot(x); title('step (with noise)');
subplot 423
[idx,score_vector,score]=nt_split(x); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on mean')
subplot 424
[idx,score_vector,score]=nt_split(x.^2); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on power')

% single channel, step in power
x=[randn(N,1); 2*randn(N,1)];  
subplot 413
plot(x); title('step in power of Gaussian noise')
subplot 427
[idx,score_vector,score]=nt_split(x); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on mean')
subplot 428
[idx,score_vector,score]=nt_split(x.^2); 
plot(score_vector); title('split on power')
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   


% 10 channels
a=nt_normcol(randn(N,5)*randn(5,10));
b=nt_normcol(randn(N,5)*randn(5,10));

% change in covariance
x=[a;b];  
figure(2); clf
subplot 311
plot(x); title('10 channel noise, change in covariance');
subplot 323
[idx,score_vector,score]=nt_split(x.^2); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on power')
subplot 324
[idx,score_vector,score]=nt_split(nt_xprod(x, 'lower')); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on cov (lower)')
subplot 325
[idx,score_vector,score]=nt_split(nt_xprod(x, 'nodiag')); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on cov (no diag)')
subplot 326
[idx,score_vector,score]=nt_split(nt_normrow(nt_xprod(x,'lower'))); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on cov (lower, normrow)')

% change in covariance, change in power
x=[a;b;2*b];  
figure(3); clf
subplot 311
plot(x); title('10 channel noise, change in covariance, then power');
subplot 323
[idx,score_vector,score]=nt_split(x.^2); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on power')
subplot 324
[idx,score_vector,score]=nt_split(nt_xprod(x, 'lower')); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on cov (lower)')
subplot 325
[idx,score_vector,score]=nt_split(nt_xprod(x, 'nodiag')); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on cov (nodiag)')
subplot 326
[idx,score_vector,score]=nt_split(nt_normrow(nt_xprod(x,'lower'))); 
disp([num2str(idx), 'score: ',num2str(score,'%.03f')]);   
plot(score_vector); title('split on cov (lower, normrow)')

