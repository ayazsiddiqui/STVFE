clear
clc
close all

cd(fileparts(mfilename('fullpath')));


%% initialize RGP
rng(1);

% altitudes
altitudes = 0:10:100;

% make class object
rgp = GP.RecursiveGaussianProcess('squaredExponential',...
    'zeroMean',altitudes);

rgp.spatialCovAmp       = 1;
rgp.spatialLengthScale  = 20;
rgp.noiseVariance       = 1e-3;

rgp.spatialCovMat = rgp.makeSpatialCovarianceMatrix(altitudes);
rgp.meanFnVector  = rgp.meanFunction(altitudes);

% guassian process
gp = GP.GaussianProcess('squaredExponential','alwaysOne','zeroMean');

gp.spatialCovAmp       = rgp.spatialCovAmp;
gp.spatialLengthScale  = rgp.spatialLengthScale;
gp.noiseVariance       = rgp.noiseVariance;

%% generate synthetic flow data

gp2 = GP.GaussianProcess('squaredExponential','exponential','windPowerLaw');
gp2.spatialCovAmp       = rgp.spatialCovAmp;
gp2.spatialLengthScale  = rgp.spatialLengthScale;
gp2.temporalCovAmp      = 1;
gp2.temporalLengthScale = 10;
gp2.noiseVariance       = rgp.noiseVariance;

% number of altitudes
nAlt = numel(altitudes);
% final time for data generation in minutes
tFinData = 120;
% time step for synthetic data generation
timeStepSynData = 1;
% standard deviation for synthetic data generation
stdDevSynData = 0.5;
% get the time series object
[synFlow,synAlt] = gp2.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData,'temporalLengthScale',1);

%% regression using traditional GP and RGP
% time step
dt = 0.25;
% algorithm final time
algFinTime = 30;
% sampling time vector
tSamp = 0:dt:algFinTime;
% number of samples
nSamp = numel(tSamp);

% preallocat sampling matrices
xSamp  = NaN(1,nSamp);
ySamp  = NaN(nSamp,1);
XTSamp = NaN(2,nSamp);

nfinAlt = 1*numel(altitudes);
finAlt = linspace(min(altitudes),max(altitudes),nfinAlt);

% preallocate matrices for GP
predMeansGP = NaN(nfinAlt,nSamp);
postVarsGP  = NaN(nfinAlt,nSamp);
stdDevGP    = NaN(nfinAlt,nSamp);
upBoundGP   = NaN(nfinAlt,nSamp);
loBoundGP   = NaN(nfinAlt,nSamp);

% preallocate matrices for KFGP
predMeansRGP = NaN(nfinAlt,nSamp);
postVarsRGP  = NaN(nfinAlt,nSamp);
stdDevRGP    = NaN(nfinAlt,nSamp);
upBoundRGP   = NaN(nfinAlt,nSamp);
loBoundRGP   = NaN(nfinAlt,nSamp);

% number of std deviations for bounds calculations
numStdDev = 1;


% make for loop
for ii = 1:nSamp
    % go to xSamp
    if ii == 1
        xSamp(ii) = altitudes(randperm(nAlt,1));
    else
        [~,maxVarIdx] = max(postVarsRGP(:,ii-1));
        xSamp(ii) = finAlt(maxVarIdx);
    end
    % measure flow at xSamp(ii) at tSamp(ii)
    fData = resample(synFlow,tSamp(ii)*60).Data;
    hData = resample(synAlt,tSamp(ii)*60).Data;
    ySamp(ii) = interp1(hData,fData,xSamp(ii));
    % augment altitude and height in XTsamp
    XTSamp(:,ii) = [xSamp(ii);tSamp(ii)];
    % recursion
    if ii == 1
        % KFGP: initial state estimate
        muGt_1   = rgp.meanFnVector;
        cGt_1    = rgp.spatialCovMat;
        % GP: covariance matrix
        covMat = gp.makeTotalCovarianceMatrix(XTSamp(:,ii));        
    else
        % KFGP: initial state estimate
        muGt_1   = predMean';
        cGt_1    = postVarMat;
        % GP: covariance matrix
        covMat = gp.augmentCovarianceMatrix(XTSamp(:,1:ii-1),...
            XTSamp(:,ii),covMat);
    end
    % KFGP: calculate prediction mean and posterior variance
    [predMean,postVarMat] =...
        rgp.calcPredMeanAndPostVar(muGt_1,cGt_1,xSamp(ii),ySamp(ii));  
    % KFGP: store them
    predMeansRGP(:,ii) = predMean;
    postVarsRGP(:,ii)  = diag(postVarMat);
    % KFGP: calculate bounds
    stdDevRGP(:,ii) = postVarsRGP(:,ii).^0.5;
    % KFGP: upper bounds = mean + x*(standard deviation)
    upBoundRGP(:,ii) = predMeansRGP(:,ii) + numStdDev*stdDevRGP(:,ii);
    % KFGP: lower bounds = mean - x*(standard deviation)
    loBoundRGP(:,ii) = predMeansRGP(:,ii) - numStdDev*stdDevRGP(:,ii);
        
    % GP: calculate prediction mean and posterior variance
    [muGP,sigGP] = ...
        gp.calcPredMeanAndPostVar(covMat,XTSamp(:,1:ii),ySamp(1:ii),...
        [finAlt;tSamp(ii)*ones(1,nfinAlt)]);
    % GP: store them
    predMeansGP(:,ii) = muGP;
    postVarsGP(:,ii)  = sigGP;
    % GP: calculate bounds
    stdDevGP(:,ii) = postVarsGP(:,ii).^0.5;
    % GP: upper bounds = mean + x*(standard deviation)
    upBoundGP(:,ii) = predMeansGP(:,ii) + numStdDev*stdDevGP(:,ii);
    % GP: lower bounds = mean - x*(standard deviation)
    loBoundGP(:,ii) = predMeansGP(:,ii) - numStdDev*stdDevGP(:,ii);
    
end

%% convert results to time series and store in strcut
regressionRes(1).predMean  = timeseries(predMeansRGP,tSamp*60);
regressionRes(1).loBound   = timeseries(loBoundRGP,tSamp*60);
regressionRes(1).upBound   = timeseries(upBoundRGP,tSamp*60);
regressionRes(1).dataSamp  = timeseries([xSamp;ySamp'],tSamp*60);
regressionRes(1).dataAlts  = timeseries(repmat(finAlt(:),1,nSamp),tSamp*60);
regressionRes(1).legend    = 'RGP';

regressionRes(2).predMean  = timeseries(predMeansGP,tSamp*60);
regressionRes(2).loBound   = timeseries(loBoundGP,tSamp*60);
regressionRes(2).upBound   = timeseries(upBoundGP,tSamp*60);
regressionRes(2).dataSamp  = timeseries([xSamp;ySamp'],tSamp*60);
regressionRes(2).dataAlts  = timeseries(repmat(finAlt(:),1,nSamp),tSamp*60);
regressionRes(2).legend    = 'GP';


%% plot the data
F = animatedPlot(synFlow,synAlt,'plotTimeStep',0.25,...
    'regressionResults',regressionRes);
