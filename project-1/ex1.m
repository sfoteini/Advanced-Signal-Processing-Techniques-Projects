%% Exercise 1: 
% Bispectrum estimation using the indirect and the direct method
%
% Foteini Savvidou

clc;
clear;
close all;

%% Question 1: 
% Construction of the process X
lambda = [0.12, 0.30, 0.42, 0.19, 0.17, 0.36];
w = 2*pi*lambda;
% low and upper limit of the uniform distribution
llimit = 0;
ulimit = 2*pi;
phi = (ulimit-llimit)*rand(1,6)+llimit;
phi(3) = phi(1) + phi(2);
phi(6) = phi(4) + phi(5);
% data length
N = 8192;
% data
k = (0:N-1)';
X = sum(cos(k*w+phi),2);
figure();
plot(X);
title('Real discrete process X');

%% Question 2: 
% Estimation of autocorrelation (128 max shiftings) and power spectrum
muX = mean(X);
maxlag = 128;
corrX = xcorr(X,maxlag);
covariance = corrX-muX^2;
% Fourier transform
fs = 1;                          % sampling frequency
fftX = fft(covariance);
k = length(fftX);
f = (0:k-1)*(fs/k);              % frequency range
powerspect = abs(fftX);
figure();
stem(f,powerspect,'.');
xlabel('Frequency');
ylabel('Magnitude');
title('Power Spectrum');

%% Question 3a: 
% Estimation of bispectrum using the indirect method with K=32, M=256
% and L=64 max shiftings for the third order cumulants
K = 32;
M = 256;
L = 64;

% Split data
subsetsX = reshape(X,[M,K]);
% Bispectrum estimation using the bispeci function from the HOSA toolbox
% i. Rectangular window
figure();
bispecInd1 = bispeci(subsetsX,L,M,0,'unbiased',128,1);
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','#D95319');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','#D95319');                % f2=0
title('Bispectrum estimated via the indirect method - Rectangular window');
legend('Bispectrum','Principal Region');

% ii. Parzen window
figure();
bispecInd2 = bispeci(subsetsX,L,M,0,'unbiased',128,0);
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','#D95319');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','#D95319');                % f2=0
title('Bispectrum estimated via the indirect method - Parzen window');
legend('Bispectrum','Principal Region');

%% Question 3b: 
% Estimation of bispectrum using the direct method with K=32, M=256 and J=0
J = 0;
D = 2*J+1;
figure();
bispecDir1 = bispecd(subsetsX,M,D,M,0);
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','#D95319');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','#D95319');                % f2=0
title('Bispectrum estimated via the direct method');
legend('Bispectrum','Principal Region');

%% Question 7a:
% Repeat steps 3a and 3b for different segment length

%% Task i. 
% Estimation of bispectrum using the indirect method with K=16, M=512
% and L=64 max shiftings for the third order cumulants 
K = 16;
M = 512;

% Split data
subsetsX = reshape(X,[M,K]);
% Rectangular window
figure();
bispecInd3 = bispeci(subsetsX,L,M,0,'unbiased',128,1);
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','#D95319');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','#D95319');                % f2=0
title('Bispectrum estimated via the indirect method - Rectangular window');
legend('Bispectrum','Principal Region');

% Parzen window
figure();
bispecInd4 = bispeci(subsetsX,L,M,0,'unbiased',128,0);
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','#D95319');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','#D95319');                % f2=0
title('Bispectrum estimated via the indirect method - Parzen window');
legend('Bispectrum','Principal Region');

% Estimation of bispectrum using the direct method with K=16, M=512 and J=0
figure();
bispecDir2 = bispecd(subsetsX,M,D,M,0);
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','#D95319');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','#D95319');                % f2=0
title('Bispectrum estimated via the direct method');
legend('Bispectrum','Principal Region');

%% Task ii. 
% Estimation of bispectrum using the indirect method with K=64, M=128
% and L=64 max shiftings for the third order cumulants 
K = 64;
M = 128;

% Split data
subsetsX = reshape(X,[M,K]);
% Rectangular window
figure();
bispecInd5 = bispeci(subsetsX,L,M,0,'unbiased',128,1);
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','#D95319');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','#D95319');                % f2=0
title('Bispectrum estimated via the indirect method - Rectangular window');
legend('Bispectrum','Principal Region');

% Parzen window
figure();
bispecInd6 = bispeci(subsetsX,L,M,0,'unbiased',128,0);
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','#D95319');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','#D95319');                % f2=0
title('Bispectrum estimated via the indirect method - Parzen window');
legend('Bispectrum','Principal Region');

% Estimation of bispectrum using the direct method with K=64, M=128 and J=0
figure();
bispecDir3 = bispecd(subsetsX,M,D,M,0);
% Display the primary area
hold on;
plot([0,0.25],[0,0.25],'Color','#D95319');            % f1=f2
plot([0.25,0.5],[0.25,0],'Color','#D95319');          % f1+f2=0.5
plot([0,0.5],[0,0],'Color','#D95319');                % f2=0
title('Bispectrum estimated via the direct method');
legend('Bispectrum','Principal Region');

%% Question 7b:
% Repeat the steps 3a and 3b for 50 realizations of the X[k] and compare
% the mean values of the power spectrum and the bispectrum
K = 32;
M = 256;
R = 50;
meanC2 = zeros(length(fftX),1);
meanC3Ind1 = zeros(M,M);
meanC3Ind2 = zeros(M,M);
meanC3Dir = zeros(M,M);

figure();
set(0,'DefaultFigureVisible','off');
for i=1:R
    % Construction of the X process
    phi = (ulimit-llimit)*rand(1,6)+llimit;
    phi(3) = phi(1) + phi(2);
    phi(6) = phi(4) + phi(5);
    % data
    k = (1:N)';
    X = sum(cos(k*w+phi),2);

    % Estimation of autocorrelation (128 max shiftings) and power spectrum
    muX = mean(X);
    corrX = xcorr(X,maxlag);
    covariance = corrX-muX^2;
    % Fourier transform
    fftX = fft(covariance);
    powerspect = abs(fftX);
    meanC2 = meanC2 + powerspect;

    % Bispectrum estimation using indirect method
    % Split data
    subsetsX = reshape(X,[M,K]);
    % i. Rectangular window
    bispecInd = bispeci(subsetsX,L,M,0,'unbiased',128,1);
    meanC3Ind1 = meanC3Ind1 + bispecInd;
    % ii. Parzen window
    bispecInd = bispeci(subsetsX,L,M,0,'unbiased',128,0);
    meanC3Ind2 = meanC3Ind2 + bispecInd;

    % Estimation of bispectrum using the direct method
    bispecDir = bispecd(subsetsX,M,D,M,0);
    meanC3Dir = meanC3Dir + bispecDir;
end

meanC2 = meanC2/R;
meanC3Ind1 = meanC3Ind1/R;
meanC3Ind2 = meanC3Ind2/R;
meanC3Dir = meanC3Dir/R;

nfft = 256;
if rem(nfft,2) == 0
    waxis = (-nfft/2:(nfft/2-1))/nfft;
else
    waxis = (-(nfft-1)/2:(nfft-1)/2)/nfft;
end

% Plot the power spectrum
fs = 1;
k = length(meanC2);
f = (0:k-1)*(fs/k);
set(0,'DefaultFigureVisible','on');
set(gcf,'Name','');
plot(f,meanC2);
xlabel('Frequency');
ylabel('Magnitude');
title('Mean value of power spectrum estimations');

% Plot the bispectrum (indirect method, rectangular window)
figure();
contour(waxis,waxis,abs(meanC3Ind1),4); grid on;
title(['Mean estimation of bispectrum estimated via the indirect ' ...
    'method - Rectangular window']);
xlabel('f1');
ylabel('f2');

% Plot the bispectrum (indirect method, Parzen window)
figure();
contour(waxis,waxis,abs(meanC3Ind2),4); grid on;
title(['Mean estimation of bispectrum estimated via the indirect ' ...
    'method - Parzen window']);
xlabel('f1');
ylabel('f2');

% Plot the bispectrum (direct method)
figure();
contour(waxis,waxis,abs(meanC3Dir),4); grid on;
title('Mean estimation of bispectrum estimated via the direct method');
xlabel('f1');
ylabel('f2');