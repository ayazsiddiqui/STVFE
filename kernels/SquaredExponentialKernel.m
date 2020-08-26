function val = SquaredExponentialKernel(p1,p2,hyperParams,covAmp)
%SQUAREDEXPONENTIALKERNEL Calculate covariance using squared exponential kernel
%   Inputs :
%       p1 - point 1
%       p2 - point 2
%       hyperParams - hyper parameters
%       covAmp - covariance amplitude

val = covAmp*exp(-0.5*sum((p1(:)-p2(:)).^2./hyperParams(:).^2));

end

