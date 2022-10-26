function [x,v] = generateX(N,q,b)
%GENERATEX Constructs a real discrete signal x[k], k = 1, ..., N, which is
% derived as the output of a MA-q process with coefficients b, driven by
% white non-Gaussian noise v[k], which is derived from an exponential
% distribution with mean value of 1.
% Inputs:
%   - N : number of samples
%   - q : the order of the MA process
%   - b : the coefficients of the MA process

    % White non-Gaussian noise (input)
    v = exprnd(1,[1,N]);
    v = v - mean(v);
    % Output of the MA-q process
    
    x = zeros(1,N);
    for k=1:N
       for i=0:q
           if k>i
               x(k) = x(k) + b(i+1)*v(k-i);
           end
       end
    end
end