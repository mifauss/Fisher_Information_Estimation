function [fy] = fY(snr,y,mode)
% Compute the pdf of Y numerically
% pdf of Y = \sqrt{\snr} X + Z is the convolution of fX and fZ
% Input:
%   snr: SNR
%   y: realization of Y, can be a scalar or a vector
%   mode: 0 - Gaussian input, 1 - binary input
% Outputs:
%   fy: fY, the pdf of Y

% Gaussian noise
fZ = @(z) 1/sqrt(2*pi)*exp(-z.^2/2);

if mode == 0
    % Gaussian input
    fX = @(x) 1/sqrt(2*pi)*exp(-x.^2/2);
    py = size(y,2);
    for cnt = 1:py
        fun = @(x) fX(x).*fZ(y(cnt) - sqrt(snr)*x);
        fy(cnt) = integral(fun,-Inf,Inf);
    end
else
    % Binary input
    p = 0.5;
    fy = fZ(y-sqrt(snr))*(1-p)+ fZ(y+sqrt(snr))*p;
end

end