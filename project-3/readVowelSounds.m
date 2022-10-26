function [fileA,fileE,fileI,fileO,fileU,Fs] = readVowelSounds(foldername)
%READVOWELSOUNDS Loads the 5 audio files of a person and keeps the first
% channel of the audio file.

    [fileA,Fs] = audioread(strcat(foldername,'A.wav'));
    fileA = fileA(:,1);
    fileE = audioread(strcat(foldername,'E.wav'));
    fileE = fileE(:,1);
    fileI = audioread(strcat(foldername,'I.wav'));
    fileI = fileI(:,1);
    fileO = audioread(strcat(foldername,'O.wav'));
    fileO = fileO(:,1);
    fileU = audioread(strcat(foldername,'U.wav'));
    fileU = fileU(:,1);
end