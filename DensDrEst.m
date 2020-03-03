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