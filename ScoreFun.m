function [f,df,rho] = ScoreFun(snr,y,mode)
% numerical score function $\rho = f'_Y/f_Y$
% Input:
%   snr: SNR
%   y: realization of Y, can be a scalar or a vector
%   mode: 0 - Gaussian input, 1 - binary input
% Outputs:
%   f: fY, the pdf of Y
%   df: dfY, the derivative of fY
%   rho: rho, the score function

% pdf f_Y and its derivative
f = fY(snr,y,mode);
df = dfY(snr,y,mode);

% score function
rho = df./f;

end