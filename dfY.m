function [dfy] = dfY(snr,y,mode)
% derivative of fY
% Input:
%   snr: SNR
%   y: realization of Y, can be a scalar or a vector
%   mode: 0 - Gaussian input, 1 - binary input
% Outputs:
%   df: dfY, the derivative of fY

% dfy = [diff(fY(snr,y),l,2)./diff(y),0];
% dfy = gradient(fY(snr,y,mode))./gradient(y);
f = @(t) fY(snr,t,mode);
dfy = (f(y+1E-8) - f(y)) / 1E-8;

end