function plotEstimations(x,x_est,N,q)
%PLOTESTIMATIONS Creates a plot of the original and the estimated signal
% x[k].
% Inputs:
%   - x    : original signal
%   - x_est: estimated signal using Giannakis' formula
%   - N    : number of samples
%   - q    : order of the MA process

    figure();
    plot(1:N,x);
    hold on;
    plot(1:N,x_est);
    title(sprintf("Original and estimated signal using Giannakis' " + ...
        "formula and q=%d",q));
    legend('Original signal','Estimated signal');
end