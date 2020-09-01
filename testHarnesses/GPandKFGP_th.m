clear
clc
close all

cd(fileparts(mfilename('fullpath')));


%% initialize KFGP
rng(1);

% altitudes
altitudes = 0:10:100;
kfgpTimeStep = 2;

% spatial kernel

kfgp = GP.KalmanFilteredGaussianProcess('squaredExponential','exponential',...
    'zeroMean',altitudes,kfgpTimeStep);

kfgp.spatialCovAmp       = 1;
kfgp.spatialLengthScale  = 20;
kfgp.temporalCovAmp      = 1;
kfgp.temporalLengthScale = 10;
kfgp.noiseVariance       = 1e-3;

kfgp.initVals = kfgp.initializeKFGP;
kfgp.spatialCovMat = kfgp.makeSpatialCovarianceMatrix(altitudes);
kfgp.spatialCovMatRoot = kfgp.calcSpatialCovMatRoot;


% guassian process
gp = GP.GaussianProcess('squaredExponential','exponential','zeroMean');

gp.spatialCovAmp       = kfgp.spatialCovAmp;
gp.spatialLengthScale  = kfgp.spatialLengthScale;
gp.temporalCovAmp      = kfgp.temporalCovAmp;
gp.temporalLengthScale = kfgp.temporalLengthScale;
gp.noiseVariance       = kfgp.noiseVariance;

%% generate synthetic flow data
% number of altitudes
nAlt = numel(altitudes);
% final time for data generation in minutes
tFinData = 120;
% time step for synthetic data generation
timeStepSynData = 1;
% standard deviation for synthetic data generation
stdDevSynData = 0.5;
% get the time series object
[synFlow,synAlt] = kfgp.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);

%% regression using traditional GP
% algorithm final time
algFinTime = 30;
% sampling time vector
tSamp = 0:kfgp.kfgpTimeStep:algFinTime;
% number of samples
nSamp = numel(tSamp);

% preallocat sampling matrices
xSamp  = NaN(1,nSamp);
ySamp  = NaN(nSamp,1);
XTSamp = NaN(2,nSamp);

% preallocate matrices for GP
predMeansGP = NaN(nAlt,nSamp);
postVarsGP  = NaN(nAlt,nSamp);
stdDevGP    = NaN(nAlt,nSamp);
upBoundGP   = NaN(nAlt,nSamp);
loBoundGP   = NaN(nAlt,nSamp);

% preallocate matrices for KFGP
predMeansKFGP = NaN(nAlt,nSamp);
postVarsKFGP  = NaN(nAlt,nSamp);
stdDevKFGP    = NaN(nAlt,nSamp);
upBoundKFGP   = NaN(nAlt,nSamp);
loBoundKFGP   = NaN(nAlt,nSamp);

% number of std deviations for bounds calculations
numStdDev = 1;


% make for loop
for ii = 1:nSamp
    % go to xSamp
    if ii == 1
        xSamp(ii) = altitudes(randperm(nAlt,1));
    else
        [~,maxVarIdx] = max(postVarsKFGP(:,ii-1));
        xSamp(ii) = altitudes(maxVarIdx);
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
        sk_k   = kfgp.initVals.s0;
        ck_k   = kfgp.initVals.sig0Mat;
        % GP: covariance matrix
        covMat = gp.makeTotalCovarianceMatrix(XTSamp(:,ii));        
    else
        % KFGP: initial state estimate
        sk_k   = skp1_kp1;
        ck_k   = ckp1_kp1;
        % GP: covariance matrix
        covMat = gp.augmentCovarianceMatrix(XTSamp(:,1:ii-1),...
            XTSamp(:,ii),covMat);
    end
    % KFGP: calculate kalman states
    [F_t,sigF_t,skp1_kp1,ckp1_kp1] = ...
        kfgp.calcKalmanStateEstimates(sk_k,ck_k,xSamp(ii),ySamp(ii));
    % KFGP: calculate prediction mean and posterior variance
    [muKFGP,sigKFGP] = kfgp.calcPredMeanAndPostVar(altitudes,F_t,sigF_t);  
    % KFGP: store them
    predMeansKFGP(:,ii) = muKFGP;
    postVarsKFGP(:,ii)  = sigKFGP;
    % KFGP: calculate bounds
    stdDevKFGP(:,ii) = postVarsKFGP(:,ii).^0.5;
    % KFGP: upper bounds = mean + x*(standard deviation)
    upBoundKFGP(:,ii) = predMeansKFGP(:,ii) + numStdDev*stdDevKFGP(:,ii);
    % KFGP: lower bounds = mean - x*(standard deviation)
    loBoundKFGP(:,ii) = predMeansKFGP(:,ii) - numStdDev*stdDevKFGP(:,ii);
        
    % GP: calculate prediction mean and posterior variance
    [muGP,sigGP] = ...
        gp.calcPredMeanAndPostVar(covMat,XTSamp(:,1:ii),ySamp(1:ii),...
        [altitudes;tSamp(ii)*ones(size(altitudes))]);
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
regressionRes(1).predMean  = timeseries(predMeansKFGP,tSamp*60);
regressionRes(1).loBound   = timeseries(loBoundKFGP,tSamp*60);
regressionRes(1).upBound   = timeseries(upBoundKFGP,tSamp*60);
regressionRes(1).dataSamp  = timeseries([xSamp;ySamp'],tSamp*60);
regressionRes(1).legend    = 'KFGP';

regressionRes(2).predMean  = timeseries(predMeansGP,tSamp*60);
regressionRes(2).loBound   = timeseries(loBoundGP,tSamp*60);
regressionRes(2).upBound   = timeseries(upBoundGP,tSamp*60);
regressionRes(2).dataSamp  = timeseries([xSamp;ySamp'],tSamp*60);
regressionRes(2).legend    = 'GP';


%% plot the data
F = animatedPlot(synFlow,synAlt,'plotTimeStep',0.25,...
    'regressionResults',regressionRes);


