function [val,varargout] = generateSyntheticFlowData(obj,altitudes,...
    finalTime,stdDev,varargin)
% parse inputs
pp = inputParser;
addParameter(pp,'timeStep',1,@isnumeric);
addParameter(pp,'meanFunc',@(x) 0);
addParameter(pp,'spatialLengthScale',obj.spatialLengthScale);
addParameter(pp,'temporalLengthScale',obj.temporalLengthScale);

parse(pp,varargin{:});

% use the gaussian process classdef to calculate the covariances
obj.spatialCovAmp = 1;
obj.spatialLengthScale = pp.Results.spatialLengthScale;
obj.temporalCovAmp = 1;
obj.temporalLengthScale = pp.Results.temporalLengthScale;

% local variables
noiseVar = 0.0001;
timeVals = 0:pp.Results.timeStep:finalTime;

% IDK why this is here, but Ben has it so I am leaving it
zstep = altitudes(2) - altitudes(1);
tstep = pp.Results.timeStep;
altitudes = (altitudes(1)-5*zstep):zstep:(altitudes(end)+5*zstep);
timeVals = (timeVals(1)-5*tstep):tstep:(timeVals(end)+5*tstep);

% altitude covariances
spatialCovMat = obj.makeSpatialCovarianceMatrix(altitudes);
spatialCovMat = spatialCovMat + noiseVar*eye(numel(altitudes));
Lz = chol(spatialCovMat);

% time covariances
temporalCovMat = obj.makeTemporalCovarianceMatrix(timeVals);
temporalCovMat = temporalCovMat + noiseVar*eye(numel(timeVals));
Lt = chol(temporalCovMat);

% random sampling
samp = stdDev*(randn(numel(altitudes),numel(timeVals)));

% mean functions
[~,M] = meshgrid(timeVals,pp.Results.meanFunc(altitudes));

% output
filterSamp = (Lz*(Lt*samp')') + M;
filterSamp = filterSamp(6:end-5,6:end-5);

% convert time to seconds and output a time series object
timeInSec = timeVals(6:end-5)*60;
altitudes = altitudes(6:end-5);
val = timeseries(filterSamp,timeInSec,'Name','SyntheticFlowData');

% extra outputs
varargout{1} = ...
    timeseries(repmat(altitudes(:),1,1,2),[timeInSec(1) timeInSec(end)]);

end

