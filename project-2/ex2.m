%% Exercise 2: 
% Validity check of Giannakis' formula
%
% Foteini Savvidou

clc;
clear;
close all;

%% Parameters of the MA process
N = 2048;       % number of samples
q = 5;          % order of the MA system
% Coefficients of the MA process
b = [1.0, 0.93, 0.85, 0.72, 0.59, -0.1];

%% Construct the signal
[x,v] = generateX(N,q,b);

%% Question 1:
% Justify the non-Gaussian character of v by calculating its skewness
mu_v = mean(v);
std_v = std(v);
skewness_v = sum((v-mu_v).^3)/((N-1)*std_v^3);
fprintf('The skewness of the noise v[k] is: %.4f\n', skewness_v);
% The skewness is not equal to zero and thus v is non-Gaussian noise.

%% Question 2:
% Estimate and plot the 3rd-order cumulants of x using the indirect method
% with K=32, M=64, L=20
K = 32;
M = 64;
L = 20;
[~,~,cum3_x,~] = bisp3cum(x,M,L,'n','u');

% Plot the 3rd order cumulants
figure();
axis=-L:L;
surf(axis,axis,cum3_x);
title('3rd order cumulants of x[k]');
xlabel('\tau_1');
ylabel('\tau_2');

figure();
contour(axis,axis,cum3_x);
title('3rd order cumulants of x[k]');
xlabel('\tau_1');
ylabel('\tau_2');

%% Question 3:
% Estimate the impulse response of the MA system using Giannakis's formula
h = GiannakisFormula(cum3_x,q,L);
fprintf("Estimation of the impulse response using Giannakis' formula" + ...
    " and order q=%d:\n",q);
disp(h);

%% Question 4:
% Estimate the impulse response of the MA system using Giannakis's formula
% considering sub-estimation and sup-estimation of the order
% Sub-estimation of the order q
q_sub = q - 2;
h_sub = GiannakisFormula(cum3_x,q_sub,L);
fprintf("Estimation of the impulse response using Giannakis' formula" + ...
    " and order q=%d:\n",q_sub);
disp(h_sub);

% Sup-estimation of the order q
q_sup = q + 3;
h_sup = GiannakisFormula(cum3_x,q_sup,L);
fprintf("Estimation of the impulse response using Giannakis' formula" + ...
    " and order q=%d:\n",q_sup);
disp(h_sup);

%% Question 5:
% Estimate the MA-q system output, plot the original and the estimated x[k]
% and find the normalized root mean square error (NRMSE)
x_est = conv(v,h,'same')';
plotEstimations(x,x_est,N,q);
nrmse = calculateNRMSE(x,x_est,N);
fprintf('NRMSE using q = %d: %.4f\n',q,nrmse);

%% Question 6:
% Repeat Step 5 for h_sub and h_sup
% Sub-estimation of the order q
x_est_sub = conv(v,h_sub,'same')';
plotEstimations(x,x_est_sub,N,q_sub);
nrmse_sub = calculateNRMSE(x,x_est_sub,N);
fprintf('NRMSE using q = %d: %.4f\n',q_sub,nrmse_sub);
% Sup-estimation of the order q
x_est_sup = conv(v,h_sup,'same')';
plotEstimations(x,x_est_sup,N,q_sup);
nrmse_sup = calculateNRMSE(x,x_est_sup,N);
fprintf('NRMSE using q = %d: %.4f\n',q_sup,nrmse_sup);

%% Question 7:
% Consider that we add a noise source of white Gaussian noise at the output
% of the system, producing a variation in the SNR
% Repeat Steps 2, 3 and 5
snr = 30:-5:-5;
n = length(snr);            % number of records
nrmse_awgn = zeros(1,n);
for i=1:n
    % Generate the noise contaminated output
    y = awgn(x,snr(i),'measured');
    % Find the 3rd order cumulants
    [~,~,cum3_y,~] = bisp3cum(y,M,L,'n','u');
    % Estimate the impulse response Giannakis's formula and q=5
    h_awgn = GiannakisFormula(cum3_y,q,L);
    % Estimate the output and find NRMSE
    y_est = conv(v,h_awgn,'same')';
    nrmse_awgn(i) = calculateNRMSE(y,y_est,N);
end
% Plot NRMSE versus SNR
figure();
plot(snr,nrmse_awgn,'-o');
xlabel('SNR [dB]');
ylabel('NRMSE');
title('NRMSE in the estimation of the output versus SNR');

%% Question 8:
% Repeat the whole process 50 times and work with the mean values of your
% results
r = 50;
NRMSE_awgn = zeros(r,n);
NRMSE = zeros(1,r);
NRMSE_sub = zeros(1,r);
NRMSE_sup = zeros(1,r);

for i=1:r
    [x_i,v_i] = generateX(N,q,b);
    % Find the 3rd order cumulants
    [~,~,cum3_xi,~] = bisp3cum(x_i,M,L,'n','u');
    % Estimate the impulse response Giannakis's formula
    h_i = GiannakisFormula(cum3_xi,q,L);
    h_i_sub = GiannakisFormula(cum3_xi,q_sub,L);
    h_i_sup = GiannakisFormula(cum3_xi,q_sup,L);
    % Estimate the output and find NRMSE
    x_est_i = conv(v_i,h_i,'same')';
    NRMSE(i) = calculateNRMSE(x_i,x_est_i,N);
    x_est_i = conv(v_i,h_i_sub,'same')';
    NRMSE_sub(i) = calculateNRMSE(x_i,x_est_i,N);
    x_est_i = conv(v_i,h_i_sup,'same')';
    NRMSE_sup(i) = calculateNRMSE(x_i,x_est_i,N);
    for j=1:n
        % Generate the noise contaminated output
        y = awgn(x_i,snr(j),'measured');
        % Find the 3rd order cumulants
        [~,~,cum3_y,~] = bisp3cum(y,M,L,'n','u');
        % Estimate the impulse response Giannakis's formula and q=5
        h_awgn = GiannakisFormula(cum3_y,q,L);
        % Estimate the output and find NRMSE
        y_est = conv(v_i,h_awgn,'same')';
        NRMSE_awgn(i,j) = calculateNRMSE(y,y_est,N);
    end
end
% Find the mean values of NRMSE
meanNRMSE = mean(NRMSE);
stdNRMSE = std(NRMSE);
meanNRMSE_sub = mean(NRMSE_sub);
stdNRMSE_sub = std(NRMSE_sub);
meanNRMSE_sup = mean(NRMSE_sup);
stdNRMSE_sup = std(NRMSE_sup);
meanNRMSE_awgn = mean(NRMSE_awgn);

% Display the results
fprintf(['\nMean value and standard deviation of NRMSE for %d ' ...
    'iterarations:\n'],r);
fprintf('\tNRMSE for q=%d: mu=%.4f std=%.4f\n',q,meanNRMSE,stdNRMSE);
fprintf('\tNRMSE for q=%d: mu=%.4f std=%.4f\n',q_sub,meanNRMSE_sub, ...
    stdNRMSE_sub);
fprintf('\tNRMSE for q=%d: mu=%.4f std=%.4f\n',q_sup,meanNRMSE_sup, ...
    stdNRMSE_sup);

% Find the 95% confidence interval of mean NRMSE with awgn using bootstrap 
% samples
alpha = 0.05;
B = 1000;       % number of bootstrap samples
bootstrap_meanNRMSE_awgn = bootstrp(B, @mean, NRMSE_awgn);
% Calculate low and upper limit for confidence interval
ci_llimit = floor((B+1)*alpha/2);
ci_ulimit = B+1-ci_llimit;
% Sort the bootstrap mean array
bootstrap_meanNRMSE_awgn_sorted = sort(bootstrap_meanNRMSE_awgn);
% Get the confidence interval
bci_meanNRMSE_awgn = zeros(n,2);
for i = 1:n
    bci_meanNRMSE_awgn(i,1) = bootstrap_meanNRMSE_awgn_sorted(ci_llimit,i);
    bci_meanNRMSE_awgn(i,2) = bootstrap_meanNRMSE_awgn_sorted(ci_ulimit,i);
end

% Plot mean NRMSE versus SNR
figure();
plot(snr,meanNRMSE_awgn,'-o');
hold on;
plot(snr,bci_meanNRMSE_awgn(:,1),'-.','Color','#A2142F');
plot(snr,bci_meanNRMSE_awgn(:,2),'-.','Color','#A2142F');
xlabel('SNR [dB]');
ylabel('NRMSE');
title('mean NRMSE in the estimation of the output versus SNR');
legend('Mean value of NRMSE','95% confidence interval');