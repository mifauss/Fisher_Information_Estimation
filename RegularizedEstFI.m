function [Jn] = RegularizedEstFI(Yi,a0,a1,kn)
% Estimator of Fisher information 
% Inputs:
%   Yi: psnr*n matrix of samples
%   a: bandwidth of the estimator
%   kn: estimator parameter
% Output:
%   Jn: estimator of Fisher information

    psnr = size(Yi,1); % n is # samples and psnr is # points of snr
    snr = 1;
    varX = 1;

    Jn = zeros(1,psnr);
    for cnt = 1:psnr
        Y = Yi(cnt,:); % vector of samples with SNR=snr(cnt)

        fn = @(t) DensEst(t,Y,a0); % Gaussian kernel density estimator
        dfn = @(t) DensDrEst(t,Y,a1); % Density derivative estimator

        rho = @(t) sqrt(3*varX)+(2*sqrt(snr)+1)*abs(t);
        fun = @(t) min(abs(dfn(t))./fn(t),rho(t)).*abs(dfn(t));
        Jn(cnt) = integral(fun,-kn,kn);
    end

end

function [fn] = DensEst(y,Yi,a)
% Gaussian kernel density estimator
% Inputs:
%   y: realization of Y, can be a scalar or a vector
%   Yi: 1*n vector of samples
%   a: bandwidth of the estimator
% Output:
%   fn: estimator of density function
    n = size(Yi,2); % n is # samples
    K = @(u) 1/sqrt(2*pi).*exp(-u.^2/2); % Gaussian kernel
    fn = 1/(n*a).*sum(K((Yi'*ones(size(y)) - ones(n,1)*y)./a),1); % kernel density estimator
end

function [dfn] = DensDrEst(y,Yi,a)
% density derivative estimator 
% Inputs:
%   y: realization of Y, can be a scalar or a vector
%   Yi: vector of samples
%   a: bandwidth of the estimator
% Output:
%   dfn: estimator of density function
    n = size(Yi,2);
    D = @(u) -u/sqrt(2*pi).*exp(-u.^2/2); % Gaussian derivative kernel
    dfn = -1/(n*a^2).*sum(D((Yi'*ones(size(y)) - ones(n,1)*y)./a),1); % density derivative estimator
end