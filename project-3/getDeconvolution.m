function [p,h] = getDeconvolution(signal,fs,llimit,ulimit,filterlimit,figtitle)
%GETDECONVOLUTION Applies a hamming window to the voice signal taking into
% account the lower and the upper limit (parameters), computes the complex
% cepstrum, lifters the cepstrum domain signal and returns the impulse
% response and the pitch signal.
% Inputs:
%   - signal : voice signal
%   - fs : sampling freuency
%   - llimit, ulimit : lower and upper limit to apply the hamming window
%   - filterlimit : limit to apply the low-pass liftering
%   - figtitle : title of the figure
% Outputs:
%   - p : pithcing signal
%   - h : impulse response
    
    % Compute the windowed signal in the time domain
    x = signal(llimit:ulimit);
    xwin = x.*hamming(ulimit-llimit+1);
    t = (0:length(xwin)-1)/fs;
    figure();
    subplot(4,1,1);
    plot(t,x);
    hold on;
    plot(t,xwin);
    xlabel('Time (s)');
    legend('Original signal','Windowed signal');
    title(figtitle);
    
    % Computer the complex cepstrum of the windowed signal
    [ccep,nd] = cceps(xwin);
    subplot(4,1,2);
    plot(t,ccep);
    xlabel('Quefrency (s)');
    title('Complex Cepstrum');

    % Liftering
    ccep_h = zeros(length(xwin),1);
    ccep_h(1:filterlimit) = ccep(1:filterlimit);
    ccep_h(end:-1:end-filterlimit-2) = ccep(end:-1:end-filterlimit-2);
    ccep_p = ccep - ccep_h;
    
    % Impluse response and pitch signal in the time domain
    p = icceps(ccep_p,nd);
    subplot(4,1,3);
    plot(t,p);
    xlabel('Time (s)');
    title('Pitch Signal');

    h = icceps(ccep_h,nd);
    subplot(4,1,4);
    plot(t,h);
    xlabel('Time (s)');
    title('Impulse response');

    sgtitle(figtitle);
end