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