function rcepstrum = plotRealCepstrum(signal,fs,figtitle)
%PLOTREALCEPSTRUM Computes and plots the real spectrum of a signal.
% Inputs:
%   - signal   : the voice signal
%   - fs       : sampling frequency
%   - figtitle : the title of the figure
% Output:
%   - rcepstrum: the real cepstrum of the signal
    
    rcepstrum = rceps(signal);
    t = (0:length(signal)-1)/fs;
    figure();
    plot(t,rcepstrum);
    xlabel('Quefrency (s)');
    title(strcat('Real Cepstrum: ',figtitle));
end