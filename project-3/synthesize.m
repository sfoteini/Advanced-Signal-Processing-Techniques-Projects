function synthesize(p,h,fs,signame)
%Synthesize Synthesizes back a voiced signal  by convoluting the impulse 
% response and the excitation signal of the original signal.
% Inputs:
%   - p : excitation signal
%   - h : impulse response
%   - fs : sampling frequency
%   - signame : name of the signal

    signal = conv(p,h);
    signal = repmat(signal,20,1);       % replicate the signal 20 times
    prompt = "\nIf you want to hear the synthesized signal, press 1\n";
    in = input(prompt,'s');
    if in == '1'
        fprintf('Now playing "%s"...',signame);
        sound(signal,fs);
    end
    prompt = "\nIf you want to see the time-domain plot of the " + ...
        "synthesized signal, press 1\n";
    in = input(prompt,'s');
    if in == '1'
        t = (0:length(signal)-1)/fs;
        figure();
        plot(t,signal);
        xlabel('Time (s)');
        title(strcat('Synthesized signal',signame));
    end
end