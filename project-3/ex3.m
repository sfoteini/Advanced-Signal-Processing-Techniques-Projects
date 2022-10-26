%% Exercise 3: 
% Cepstrum via Homomorphic Filtering
%
% Foteini Savvidou

clc;
clear;
close all;

%% Task 1:
% Load 10 voice samples of a male and female individual making the five 
% vowel sounds (a, e, i, o, u).
female_folder_name = 'Samples\Female\';
[femaleA,femaleE,femaleI,femaleO,femaleU,fs] = readVowelSounds( ...
    female_folder_name);
male_folder_name = 'Samples\Male\';
[maleA,maleE,maleI,maleO,maleU,~] = readVowelSounds(male_folder_name);

%% Plot the signals in the time domain
prompt = "If you want to see the plots of the 10 signals in the " + ...
    "time domain, press 1\n";
in = input(prompt,'s');
if in == '1'
    plotSignal(femaleA,fs,'Vowel A - Female');
    plotSignal(femaleE,fs,'Vowel E - Female');
    plotSignal(femaleI,fs,'Vowel I - Female');
    plotSignal(femaleO,fs,'Vowel O - Female');
    plotSignal(femaleU,fs,'Vowel U - Female');
    plotSignal(maleA,fs,'Vowel A - Male');
    plotSignal(maleE,fs,'Vowel E - Male');
    plotSignal(maleI,fs,'Vowel I - Male');
    plotSignal(maleO,fs,'Vowel O - Male');
    plotSignal(maleU,fs,'Vowel U - Male');
end

%% Task 2:
% Compute the cepstrum and plot the results
rcep_f_a = plotRealCepstrum(femaleA,fs,'Vowel A - Female');
rcep_f_e = plotRealCepstrum(femaleE,fs,'Vowel E - Female');
rcep_f_i = plotRealCepstrum(femaleI,fs,'Vowel I - Female');
rcep_f_o = plotRealCepstrum(femaleO,fs,'Vowel O - Female');
rcep_f_u = plotRealCepstrum(femaleU,fs,'Vowel U - Female');
rcep_m_a = plotRealCepstrum(maleA,fs,'Vowel A - Male');
rcep_m_e = plotRealCepstrum(maleE,fs,'Vowel E - Male');
rcep_m_i = plotRealCepstrum(maleI,fs,'Vowel I - Male');
rcep_m_o = plotRealCepstrum(maleO,fs,'Vowel O - Male');
rcep_m_u = plotRealCepstrum(maleU,fs,'Vowel U - Male');

%% Task 3:
% Lifter the cepstrum domain signals.
% Design a Hamming window to remove the transfer function dependency and 
% compute the time domain signal of the corresponding windowed result to 
% obtain the deconvolved signal.

[p_f_a,h_f_a] = getDeconvolution(femaleA,fs,10001,10900,10,'Vowel A - Female');
[p_f_e,h_f_e] = getDeconvolution(femaleE,fs,10001,10950,12,'Vowel E - Female');
[p_f_i,h_f_i] = getDeconvolution(femaleI,fs,10001,10920,10,'Vowel I - Female');
[p_f_o,h_f_o] = getDeconvolution(femaleO,fs,10001,10900,15,'Vowel O - Female');
[p_f_u,h_f_u] = getDeconvolution(femaleU,fs,10001,10940,8,'Vowel U - Female');
[p_m_a,h_m_a] = getDeconvolution(maleA,fs,10001,10970,35,'Vowel A - Male');
[p_m_e,h_m_e] = getDeconvolution(maleE,fs,10001,10970,65,'Vowel E - Male');
[p_m_i,h_m_i] = getDeconvolution(maleI,fs,10001,10970,70,'Vowel I - Male');
[p_m_o,h_m_o] = getDeconvolution(maleO,fs,10001,10970,40,'Vowel O - Male');
[p_m_u,h_m_u] = getDeconvolution(maleU,fs,10001,10970,80,'Vowel U - Male');

%% Task 4:
% Synthesize back the voiced signals
synthesize(p_f_a,h_f_a,fs,'Vowel A - Female');
synthesize(p_f_e,h_f_e,fs,'Vowel E - Female');
synthesize(p_f_i,h_f_i,fs,'Vowel I - Female');
synthesize(p_f_o,h_f_o,fs,'Vowel O - Female');
synthesize(p_f_u,h_f_u,fs,'Vowel U - Female');
synthesize(p_m_a,h_m_a,fs,'Vowel A - Male');
synthesize(p_m_e,h_m_e,fs,'Vowel E - Male');
synthesize(p_m_i,h_m_i,fs,'Vowel I - Male');
synthesize(p_m_o,h_m_o,fs,'Vowel O - Male');
synthesize(p_m_u,h_m_u,fs,'Vowel U - Male');