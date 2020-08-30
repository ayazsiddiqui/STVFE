clear
clc
close all

cd(fileparts(mfilename('fullpath')));


%% initialize GP
rng(1);

gp = GP.GaussianProcess(1,'squaredExponential','exponential');

gp.spatialCovAmp       = 1;
gp.spatialLengthScale  = 20;
gp.temporalCovAmp      = 1;
gp.temporalLengthScale = 10;
gp.noiseVariance       = 1e-3;

%% generate synthetic flow data
% altitudes
altitudes = 0:10:100;
nAlt = numel(altitudes);
% final time for data generation in minutes
tFinData = 30;
% time step for synthetic data generation
timeStepSynData = 1;
% standard deviation for synthetic data generation
stdDevSynData = 0.5;
% get the time series object
[synFlow,synAlt] = gp.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);


%% regression using traditional GP
% sampling time step
dt = .5;
% algorithm final time
algFinTime = 30;
% sampling time vector
tSamp = 0:dt:algFinTime;
% number of samples
nSamp = numel(tSamp);

% preallocate matrices
predMeans = NaN(nAlt,nSamp);
postVars  = NaN(nAlt,nSamp);
xSamp     = NaN(1,nSamp);
ySamp     = NaN(nSamp,1);
XTSamp    = NaN(2,nSamp);

stdDev    = NaN(nAlt,nSamp);
upBound   = NaN(nAlt,nSamp);
loBound   = NaN(nAlt,nSamp);

% number of std deviations for bounds calculations
numStdDev = 1;

% make for loop
for ii = 1:nSamp
    % go to xSamp
    if ii == 1
        xSamp(ii) = altitudes(randperm(nAlt,1));
    else
        [~,maxVarIdx] = max(postVars(:,ii-1));
        xSamp(ii) = altitudes(maxVarIdx);
    end
    % measure flow at xSamp(ii) at tSamp(ii)
    fData = resample(synFlow,tSamp(ii)*60).Data;
    hData = resample(synAlt,tSamp(ii)*60).Data;
    ySamp(ii) = interp1(hData,fData,xSamp(ii));
    % augment altitude and height in XTsamp
    XTSamp(:,ii) = [xSamp(ii);tSamp(ii)];
    % add new point to GP
    if ii == 1
        covMat = gp.makeTotalCovarianceMatrix(XTSamp(:,ii));
    else
        covMat = gp.augmentCovarianceMatrix(XTSamp(:,1:ii-1),...
            XTSamp(:,ii),covMat);
    end
    % calculate prediction mean and posterior variance
    [mu,sig] = gp.calcPredMeanAndPostVar(...
        covMat,XTSamp(:,1:ii),ySamp(1:ii),...
        [altitudes;tSamp(ii)*ones(size(altitudes))]);
    % store them
    predMeans(:,ii) = mu;
    postVars(:,ii)   = sig;
    % calculate bounds
    stdDev(:,ii) = postVars(:,ii).^0.5;
    % % % upper bounds = mean + x*(standard deviation)
    upBound(:,ii) = predMeans(:,ii) + numStdDev*stdDev(:,ii);
    % % % lower bounds = mean + x*(standard deviation)
    loBound(:,ii) = predMeans(:,ii) - numStdDev*stdDev(:,ii);
    
end

% put convert results to time series and store in strcut
regressionRes(1).predMean  = timeseries(predMeans,tSamp*60);
regressionRes(1).loBound   = timeseries(loBound,tSamp*60);
regressionRes(1).upBound   = timeseries(upBound,tSamp*60);
regressionRes(1).dataSamp  = timeseries([xSamp;ySamp'],tSamp*60);


%% plot the data
F = animatedPlot(synFlow,synAlt,'regressionResults',regressionRes);




