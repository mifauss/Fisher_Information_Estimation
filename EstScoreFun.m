function [rhon,fn,dfn] = EstScoreFun(Yi,y,a0,a1)
% Estimator of Fisher information 
% Inputs:
%   Yi: psnr*n matrix of samples
%   y: realization of Y, can be a scalar or a vector
%   a: bandwidth of the estimator
% Output:
%   fn: estimator of fY, the pdf of Y
%   dfn: estimator of dfY, the derivative of fY
%   rhon: estimator of rho, the score function

    
% Gaussian kernel density estimator
fn = DensEst(y,Yi,a0);

% Density derivative estimator
dfn = DensDrEst(y,Yi,a1);

% score function estimator
rhon = dfn./(fn+(fn==0)*eps);


end