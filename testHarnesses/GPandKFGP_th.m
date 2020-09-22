clear
clc
close all

cd(fileparts(mfilename('fullpath')));


%% initialize KFGP
rng(1);

% altitudes
altitudes = 0:100:1000;
kfgpTimeStep = 0.1;

% make class object
kfgp = GP.KalmanFilteredGaussianProcess('squaredExponential','squaredExponential',...
    'zeroMean',altitudes,kfgpTimeStep);

kfgp.spatialCovAmp       = 1;
kfgp.spatialLengthScale  = 220;
kfgp.temporalCovAmp      = 1;
kfgp.temporalLengthScale = 22;
kfgp.noiseVariance       = 0.1;

kfgp.initVals = kfgp.initializeKFGP;
kfgp.spatialCovMat = kfgp.makeSpatialCovarianceMatrix(altitudes);
kfgp.spatialCovMatRoot = kfgp.calcSpatialCovMatRoot;


% guassian process
gp = GP.GaussianProcess('squaredExponential','squaredExponential','zeroMean');

gp.spatialCovAmp       = kfgp.spatialCovAmp;
gp.spatialLengthScale  = kfgp.spatialLengthScale;
gp.temporalCovAmp      = kfgp.temporalCovAmp;
gp.temporalLengthScale = kfgp.temporalLengthScale;
gp.noiseVariance       = kfgp.noiseVariance;

%% generate synthetic flow data
gp2 = GP.GaussianProcess('squaredExponential','squaredExponential','windPowerLaw');
gp2.spatialCovAmp       = kfgp.spatialCovAmp;
gp2.spatialLengthScale  = kfgp.spatialLengthScale;
gp2.temporalCovAmp      = kfgp.temporalCovAmp;
gp2.temporalLengthScale = kfgp.temporalLengthScale;
gp2.noiseVariance       = kfgp.noiseVariance;
% number of altitudes
nAlt = numel(altitudes);
% final time for data generation in minutes
tFinData = 120;
% time step for synthetic data generation
timeStepSynData = 1;
% standard deviation for synthetic data generation
stdDevSynData = 1;
% get the time series object
[synFlow,synAlt] = gp2.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);

%% regression using traditional GP
% algorithm final time
algFinTime = 30;
% sampling time vector
tSamp = 0:kfgp.kfgpTimeStep:algFinTime;
% number of samples
nSamp = numel(tSamp);

% preallocat sampling matrices
xSamp   = NaN(1,nSamp);
flowVal = NaN(1,nSamp);
ySamp   = NaN(nSamp,1);
XTSamp  = NaN(2,nSamp);

nfinAlt = 1*numel(altitudes);
finAlt = linspace(min(altitudes),max(altitudes),nfinAlt);

meanFnVec = gp2.meanFunction(finAlt);

% preallocate matrices for GP
predMeansGP = NaN(nfinAlt,nSamp);
postVarsGP  = NaN(nfinAlt,nSamp);
RMSEFitKFGP = NaN(1,nSamp);
stdDevGP    = NaN(nfinAlt,nSamp);
upBoundGP   = NaN(nfinAlt,nSamp);
loBoundGP   = NaN(nfinAlt,nSamp);

% preallocate matrices for KFGP
predMeansKFGP = NaN(nfinAlt,nSamp);
postVarsKFGP  = NaN(nfinAlt,nSamp);
RMSEFitGP     = NaN(1,nSamp);
stdDevKFGP    = NaN(nfinAlt,nSamp);
upBoundKFGP   = NaN(nfinAlt,nSamp);
loBoundKFGP   = NaN(nfinAlt,nSamp);

% number of std deviations for bounds calculations
numStdDev = 1;


% make for loop
for ii = 1:nSamp
    % go to xSamp
    if ii == 1
        xSamp(ii) = altitudes(randperm(nAlt,1));
    else
        [~,maxVarIdx] = max(postVarsKFGP(:,ii-1));
        xSamp(ii) = finAlt(maxVarIdx);
    end
    % measure flow at xSamp(ii) at tSamp(ii)
    fData = resample(synFlow,tSamp(ii)*60).Data;
    hData = resample(synAlt,tSamp(ii)*60).Data;
    flowVal(ii) = interp1(hData,fData,xSamp(ii));
    ySamp(ii) =  gp2.meanFunction(xSamp(ii)) - flowVal(ii);
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
    [muKFGP,sigKFGP] = kfgp.calcPredMeanAndPostVar(finAlt,F_t,sigF_t);  
    % KFGP: store them
    predMeansKFGP(:,ii) = meanFnVec(:) - muKFGP;
    postVarsKFGP(:,ii)  = sigKFGP;
    % percentage fit to data
    RMSEFitKFGP(ii) = mean(((predMeansKFGP(:,ii) - ...
        interp1(hData,fData,nfinAlt))/mean(meanFnVec)).^2)^0.5; 
    % KFGP: calculate bounds
    stdDevKFGP(:,ii) = postVarsKFGP(:,ii).^0.5;
    % KFGP: upper bounds = mean + x*(standard deviation)
    upBoundKFGP(:,ii) = predMeansKFGP(:,ii) + numStdDev*stdDevKFGP(:,ii);
    % KFGP: lower bounds = mean - x*(standard deviation)
    loBoundKFGP(:,ii) = predMeansKFGP(:,ii) - numStdDev*stdDevKFGP(:,ii);
        
    % GP: calculate prediction mean and posterior variance
    [muGP,sigGP] = ...
        gp.calcPredMeanAndPostVar(covMat,XTSamp(:,1:ii),ySamp(1:ii),...
        [finAlt;tSamp(ii)*ones(1,nfinAlt)]);
    % GP: store them
    predMeansGP(:,ii) = meanFnVec(:) - muGP;
    postVarsGP(:,ii)  = sigGP;
    % percentage fit to data
    RMSEFitGP(ii) = mean(((predMeansGP(:,ii) - ...
        interp1(hData,fData,nfinAlt))/mean(meanFnVec)).^2)^0.5; 
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
regressionRes(1).dataSamp  = timeseries([xSamp;flowVal],tSamp*60);
regressionRes(1).dataAlts  = timeseries(repmat(finAlt(:),1,nSamp),tSamp*60);
regressionRes(1).legend    = 'KFGP';

regressionRes(2).predMean  = timeseries(predMeansGP,tSamp*60);
regressionRes(2).loBound   = timeseries(loBoundGP,tSamp*60);
regressionRes(2).upBound   = timeseries(upBoundGP,tSamp*60);
regressionRes(2).dataSamp  = timeseries([xSamp;flowVal],tSamp*60);
regressionRes(2).dataAlts  = timeseries(repmat(finAlt(:),1,nSamp),tSamp*60);
regressionRes(2).legend    = 'GP';


%% plot the data
F = animatedPlot(synFlow,synAlt,'plotTimeStep',0.25,...
    'regressionResults',regressionRes,'wait',true);


