function [Jn] = EstFisherInfo(Yi,a0,a1,kn)
% Estimator of Fisher information 
% Inputs:
%   Yi: psnr*n matrix of samples
%   a: bandwidth of the estimator
%   kn: estimator parameter
% Output:
%   Jn: estimator of Fisher information

psnr = size(Yi,1); % n is # samples and psnr is # points of snr

Jn = zeros(1,psnr);
for cnt = 1:psnr
    Y = Yi(cnt,:); % vector of samples with SNR=snr(cnt)
    
    fn = @(t) DensEst(t,Y,a0); % Gaussian kernel density estimator
    dfn = @(t) DensDrEst(t,Y,a1); % Density derivative estimator
    fun = @(t) (dfn(t)).^2./(fn(t)+(fn(t)==0)*eps);
    Jn(cnt) = integral(fun,-kn,kn);
        
%     rho = EstScoreFun(Y,Y,a);
%     J2(cnt) = mean(rho.^2);
%     J3(cnt) = mean( min(rho.^2, (sqrt(3)+(2*sqrt(snr)+1).*Y).^2) );
end



end