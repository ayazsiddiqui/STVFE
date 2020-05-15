function [windSpeedOut,out] = genWindv2(z,lz,t,lt,std,varargin)
% Function generates wind data with given length scales and standard
% deviation h around a given mean function m(z)
% Inputs
%     z - vector of all altitude values
%     clz - characteristic length in z
%     t - vector of all time values
%     clt - characteristic length in t
%     h - standard deviation around mean
%     m - mean function evaluated at all z values
% Outputs
%     windSpeed - wind speed (each time step is a column)
%     Z - z values in a grid
%     T - t values in a grid
    

zstep = z(2)-z(1);
tstep = t(2)-t(1);
z = (z(1)-5*zstep):zstep:(z(end)+5*zstep);
t = (t(1)-5*tstep):tstep:(t(end)+5*tstep);
z = z(:);
t = t(:);

covz = exp(-(z-z').^2/(2*lz^2));
epz = .0001*diag(ones(1,length(z)));
covz = covz + epz;
Lz = chol(covz);

covt = exp(-(t-t').^2/(2*lt^2));
ept = .0001*diag(ones(1,length(t)));
covt = covt + ept;
Lt = chol(covt);

if nargin == 5
    samp = std*(randn(length(z),length(t)));
    mf = @(x) 0;
elseif nargin == 6
    mf = varargin{1};
    samp = std*(randn(length(z),length(t)));
elseif nargin == 7
    mf = varargin{1};
    samp = varargin{2};
end

[~,M] = meshgrid(t,mf(z));

filterSamp = (Lz*(Lt*samp')') + M;
windSpeedOut = filterSamp(6:end-5,6:end-5);
out = windSpeedOut;

end