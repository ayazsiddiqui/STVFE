clear
clc
close all

cd(fileparts(mfilename('fullpath')));


%% initialize RGP
rng(1);

% altitudes
altitudes = 0:10:100;
% number of altitudes
nAlt = numel(altitudes);

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

% final time for data generation in minutes
SD.tFinData = 120;
% time step for synthetic data generation
SD.timeStepSynData = 1;
% standard deviation for synthetic data generation
SD.stdDevSynData = 0.5;
% get the time series object
[synFlow,synAlt] = gp2.generateSyntheticFlowData(altitudes,SD.tFinData,SD.stdDevSynData,...
    'timeStep',SD.timeStepSynData,'temporalLengthScale',1);

%% regression using traditional GP and RGP
% time step
A.dt = 0.25;
% algorithm final time
A.algFinTime = 30;

% sampling time vector
S.tSamp = 0:A.dt:A.algFinTime;
% number of samples
nSamp = numel(S.tSamp);
% preallocat sampling matrices
S.xSamp  = NaN(1,nSamp);
S.ySamp  = NaN(nSamp,1);
S.XTSamp = NaN(2,nSamp);

% preallocate matrices for GP
GP.predMean = NaN(nAlt,nSamp);
GP.postVars = NaN(nAlt,nSamp);
GP.stdDev   = NaN(nAlt,nSamp);
GP.upBound  = NaN(nAlt,nSamp);
GP.loBound  = NaN(nAlt,nSamp);

% preallocate matrices for RGP
RGP.predMean = NaN(nAlt,nSamp);
RGP.postVars = NaN(nAlt,nSamp);
RGP.stdDev   = NaN(nAlt,nSamp);
RGP.upBound  = NaN(nAlt,nSamp);
RGP.loBound  = NaN(nAlt,nSamp);

% preallocate matrices for colored RGP
RGPC.predMean = NaN(nAlt,nSamp);
RGPC.postVars = NaN(nAlt,nSamp);
RGPC.stdDev   = NaN(nAlt,nSamp);
RGPC.upBound  = NaN(nAlt,nSamp);
RGPC.loBound  = NaN(nAlt,nSamp);

rgp.filterTimeConstant = 0.01;
rgp.sampleTime         = A.dt;
% number of std deviations for bounds calculations
A.numStdDev = 1;


% make for loop
for ii = 1:nSamp
    % go to xSamp
    if ii == 1
        S.xSamp(ii) = altitudes(randperm(nAlt,1));
    else
%         [~,maxVarIdx] = max(RGP.postVars(:,ii-1));
        S.xSamp(ii) = altitudes(randperm(nAlt,1));
    end
    % measure flow at xSamp(ii) at tSamp(ii)
    S.ySamp(ii) = interp1(resample(synAlt,S.tSamp(ii)*60).Data,...
        resample(synFlow,S.tSamp(ii)*60).Data,S.xSamp(ii));
    % augment altitude and height in XTsamp
    S.XTSamp(:,ii) = [S.xSamp(ii);S.tSamp(ii)];
    % recursion
    if ii == 1
        % RGP: initial state estimate
        muGt_1   = rgp.meanFnVector;
        cGt_1    = rgp.spatialCovMat;
        % GP: covariance matrix
        covMat = gp.makeTotalCovarianceMatrix(S.XTSamp(:,ii));     
        % colored RGP: 
        CmuGt_1   = rgp.meanFnVector;
        CcGt_1    = rgp.spatialCovMat;
    else
        % RGP: initial state estimate
        muGt_1   = predMean';
        cGt_1    = postVarMat;
        % GP: covariance matrix
        covMat = gp.augmentCovarianceMatrix(S.XTSamp(:,1:ii-1),...
            S.XTSamp(:,ii),covMat);
        % colored RGP: initial state estimate
        CmuGt_1   = CpredMean';
        CcGt_1    = CpostVarMat;
    end
    % RGP: calculate prediction mean and posterior variance
    [predMean,postVarMat] =...
        rgp.calcPredMeanAndPostVar(muGt_1,cGt_1,S.xSamp(ii),S.ySamp(ii));  
    % RGP: store them
    RGP.predMean(:,ii) = predMean;
    RGP.postVars(:,ii)  = diag(postVarMat);
    % RGP: calculate bounds
    RGP.stdDev(:,ii) = RGP.postVars(:,ii).^0.5;
    % RGP: upper bounds = mean + x*(standard deviation)
    RGP.upBound(:,ii) = RGP.predMean(:,ii) + A.numStdDev*RGP.stdDev(:,ii);
    % RGP: lower bounds = mean - x*(standard deviation)
    RGP.loBound(:,ii) = RGP.predMean(:,ii) - A.numStdDev*RGP.stdDev(:,ii);
        
    % GP: calculate prediction mean and posterior variance
    [muGP,sigGP] = ...
        gp.calcPredMeanAndPostVar(covMat,S.XTSamp(:,1:ii),S.ySamp(1:ii),...
        [altitudes;S.tSamp(ii)*ones(1,nAlt)]);
    % GP: store them
    GP.predMean(:,ii) = muGP;
    GP.postVars(:,ii)  = sigGP;
    % GP: calculate bounds
    GP.stdDev(:,ii) = GP.postVars(:,ii).^0.5;
    % GP: upper bounds = mean + x*(standard deviation)
    GP.upBound(:,ii) = GP.predMean(:,ii) + A.numStdDev*GP.stdDev(:,ii);
    % GP: lower bounds = mean - x*(standard deviation)
    GP.loBound(:,ii) = GP.predMean(:,ii) - A.numStdDev*GP.stdDev(:,ii);
    
    % colored RGP: calculate prediction mean and posterior variance
    [CpredMean,CpostVarMat] =...
        rgp.calcColoredPredMeanAndPostVar(CmuGt_1,CcGt_1,S.xSamp(ii),S.ySamp(ii));  
    % colored RGP: store them
    RGPC.predMean(:,ii) = CpredMean;
    RGPC.postVars(:,ii)  = diag(CpostVarMat);
    % colored RGP: calculate bounds
    RGPC.stdDev(:,ii) = RGPC.postVars(:,ii).^0.5;
    % colored RGP: upper bounds = mean + x*(standard deviation)
    RGPC.upBound(:,ii) = RGPC.predMean(:,ii) + A.numStdDev*RGPC.stdDev(:,ii);
    % colored RGP: lower bounds = mean - x*(standard deviation)
    RGPC.loBound(:,ii) = RGPC.predMean(:,ii) - A.numStdDev*RGPC.stdDev(:,ii);
    
end

%% convert results to time series and store in strcut
regressionRes(1).predMean  = timeseries(RGP.predMean,S.tSamp*60);
regressionRes(1).loBound   = timeseries(RGP.loBound,S.tSamp*60);
regressionRes(1).upBound   = timeseries(RGP.upBound,S.tSamp*60);
regressionRes(1).dataSamp  = timeseries([S.xSamp;S.ySamp'],S.tSamp*60);
regressionRes(1).dataAlts  = timeseries(repmat(altitudes(:),1,nSamp),S.tSamp*60);
regressionRes(1).legend    = 'RGP';

regressionRes(2).predMean  = timeseries(GP.predMean,S.tSamp*60);
regressionRes(2).loBound   = timeseries(GP.loBound,S.tSamp*60);
regressionRes(2).upBound   = timeseries(GP.upBound,S.tSamp*60);
regressionRes(2).dataSamp  = timeseries([S.xSamp;S.ySamp'],S.tSamp*60);
regressionRes(2).dataAlts  = timeseries(repmat(altitudes(:),1,nSamp),S.tSamp*60);
regressionRes(2).legend    = 'GP';

regressionRes(3).predMean  = timeseries(RGPC.predMean,S.tSamp*60);
regressionRes(3).loBound   = timeseries(RGPC.loBound,S.tSamp*60);
regressionRes(3).upBound   = timeseries(RGPC.upBound,S.tSamp*60);
regressionRes(3).dataSamp  = timeseries([S.xSamp;S.ySamp'],S.tSamp*60);
regressionRes(3).dataAlts  = timeseries(repmat(altitudes(:),1,nSamp),S.tSamp*60);
regressionRes(3).legend    = 'CRGP';

%% plot the data
t = tiledlayout('flow');
nexttile;
semilogy(1:nSamp,vecnorm(GP.predMean-RGP.predMean));
hold on
semilogy(1:nSamp,vecnorm(GP.predMean-RGPC.predMean),'-.o');
semilogy(1:nSamp,vecnorm(RGP.predMean-RGPC.predMean),'-+');
legend('GP vs. RPG','GP vs. RGPC','RGP vs. RGPC');
ylabel('Pred. mean difference')


%% amimate the data
F = animatedPlot(synFlow,synAlt,'plotTimeStep',0.25,...
    'regressionResults',regressionRes);
