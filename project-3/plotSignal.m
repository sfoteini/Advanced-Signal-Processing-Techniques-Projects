function plotSignal(signal,fs,figtitle)
%PLOTSIGNAL Plots the given signal in the time domain.
% Inputs:
%   - signal : a voice signal
%   - fs     : sampling frequency
%   - title  : the title of the generated figure

    t = (0:length(signal)-1)/fs;
    figure();
    plot(t,signal);
    xlabel('Time (s)');
    title(figtitle);
end