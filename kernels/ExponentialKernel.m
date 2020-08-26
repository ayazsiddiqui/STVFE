function val = ExponentialKernel(p1,p2,hyperParams,covAmp)
%EXPONENTIALKERNEL Calculate covariance using exponential kernel
%   Inputs :
%       p1 - point 1
%       p2 - point 2
%       hyperParams - hyper parameters
%       covAmp - covariance amplitude

val = covAmp*exp(-sum(abs(p1(:)-p2(:))./hyperParams(:)));

end

