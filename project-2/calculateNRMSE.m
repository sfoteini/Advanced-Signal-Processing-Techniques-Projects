function nrmse = calculateNRMSE(x,x_est,N)
%CALCULATENRMSE Calculates the normalized root mean square error of two 
% signals.
% Inputs:
%   - x, x_est : two signals
%   - N        : number of samples
% Outputs:
%   - nrmse    : normalized root mean square error

    rmse = sqrt(sum((x_est-x').^2)/N);
    nrmse = rmse/(max(x)-min(x));
end