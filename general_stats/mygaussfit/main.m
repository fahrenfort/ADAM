function main
%% clearing
fclose all; close all; clear all; 
clc; pause(0.1);

%% data
x=1:100;
sigma=30; mu=40; A=3;
y=A*exp(-(x-mu).^2/(2*sigma^2))+rand(size(x))*0.3;
plot(x,y,'.');

%% fitting
[sigmaNew,muNew,Anew]=mygaussfit(x,y);
y=Anew*exp(-(x-muNew).^2/(2*sigmaNew^2));
hold on; plot(x,y,'.r');


%% fitting basis set
new_new_core = mean(subj_CTFs); % get mirrored basis set
new_new_core = new_new_core-min(new_new_core); % set min value to 0
new_new_core = padarray([new_new_core(end) new_new_core]',3)'; % add some zeros on both sides
new_new_core = new_new_core/sum(new_new_core); % normalize gaussian

y = new_new_core;
x = 1:numel(y);
figure;
plot(x,y);
[sigma,mu,A]=mygaussfit(x,y);
y = mygausswin(x,sigma,mu,A);

hold on; plot(x,y,'r');
%%
figure;
y1 = new_new_core(5:12);
x1 = x(5:12);
plot(x1,y1);
y2 = mygausswin(x1,sigma,mu,A);
hold on;
plot(x1,y2);

%% create nice gauss USE THIS FOR FEM IF NO MODEL IS SPECIFIED
nCond = 10;
%set1 = mygausswin(1:nCond,1.5,0,1);
% if mod(nCond,2) == 0
%     set1 = mygausswin(1:nCond,1.5)'; % use gaussian tuning curve
% else
%     set1 = mygausswin(1:nCond-1,1.5)';
%     set1 = [set1(end); set1];
% end
set1 = mygausswin(1:nCond,nCond/6,round(nCond/2))';
%set1 = set1/sum(set1);
%y = y/sum(y); % normalize gaussian
%x = 1:numel(y);

if mod(nCond,2) == 0
    set2 = [gausswin(nCond-1); 0]; % use gaussian tuning curve
else
    set2 = [0; gausswin(nCond-2); 0 ];
end

figure; plot(set1); hold on; plot(set2);