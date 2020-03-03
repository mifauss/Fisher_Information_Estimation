function [J] = FisherInfo(snr,mode)
% compute Fisher information
% Input:
%   snr: SNR
%   mode: 0 - Gaussian input, 1 - binary input
% Outputs:
%   J: Fisher information

psnr = size(snr,2);

J = zeros(1,psnr);

for cnt = 1:psnr
    snr0 = snr(cnt);
    
    % Fisher information
    fy = @(t)fY(snr0,t,mode);
    dfy = @(t)dfY(snr0,t,mode);
    func = @(t) (dfy(t)).^2./(fy(t)+(fy(t)==0).*eps);
    J(cnt) = integral(func,-Inf,Inf,'RelTol',1e-2,'AbsTol',1e-3);

end

end