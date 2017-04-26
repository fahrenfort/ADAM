% demo on how to obtain patterns from classifier
% weights (filters)
% Stefan Haufe, 2014

M = 20; %number of features
N = 500; %number of trials
Nplus = N/2; %number of trials in positive class


% generate data for two classes according to LDA model

% binary labels
labels = [ones(Nplus, 1); -ones(N-Nplus, 1)];

% data
% create a random covariance matrix
O = orth(randn(M));
cov_theo = O*diag(exp(2*randn(M, 1)))*O'; 
cov_theo = cov_theo ./ norm(cov_theo);
% create offset: class means differ only in first five features, not in the
% last five. This is the vector we want to recover
offset = zeros(1, M); offset(1:5) = 1*randn(1, 5); 
% generate data finally
trials = randn(N, M)*sqrtm(cov_theo) + [repmat(offset, Nplus, 1); zeros(N-Nplus, M)];
% trials = trials - repmat(mean(trials), N, 1);

% plot trials of two classes in 1st and 6th dimension
% distributions differ only in the 1st, but not in the 6th
plot(trials(1:Nplus, 1), trials(1:Nplus, 6), 'b.')
hold on
plot(trials(Nplus+1:N, 1), trials(Nplus+1:N, 6), 'r.')
grid on

% train a classifier based on least-squares regression to the class labels
% (replace this by SVM, logistic regression, etc.)
wb = [trials ones(N, 1)]\labels;
w = wb(1:M); % weight vector (filter)
b = wb(end); % offset

% obtain pattern 
a = cov(trials)*w;
%a = corrcoef(trials)*w;

% alternative way to obtain patterns, which does not explicitly compute or
% store the covariance matrix cov(trials), and is therefore more space and
% time efficient especially if M >> N. Precisely, it does not require
% more space than the actual data in trials. This is the recommended 
% way to compute patterns for fMRI data.
trials_ = (trials - repmat(mean(trials), N, 1))./sqrt(N-1);
a2 = trials_'*(trials_*w); 

% compare filters, patterns, and the original classmean differences we want
% to recover
figure;
plot(offset./std(offset), 'b')
hold on
plot(w./std(w), 'r')
plot(a./std(a), 'g')
plot(a2./std(a2), 'k--')
grid on
legend('true pattern', 'filter', 'pattern', 'pattern2')


