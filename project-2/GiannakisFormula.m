function h = GiannakisFormula(cum3,q,L)
%GIANNAKIS FORMULA Estimates the impulse response of the MA system using 
% Giannakis's formula.
% Inputs:
%   - cum3: 3rd order cumulants of the input signal
%   - q   : order of the MA process
%   - L   : max-shiftings used in the estimation of the 3rd order cumulants 
% Outputs:
%   - h   : impulse response of the MA system

    h = cum3(q+L+1,(L+1):(q+L+1))/cum3(q+L+1,L+1);
end