clear
clc
close all

cd(fileparts(mfilename('fullpath')));


%% initialize KFGP
rng(1);

% altitudes
altitudes = 0:10:100;
kfgpTimeStep = 0.05;

% spatial kernel

kfgp = GP.KalmanFilteredGaussianProcess('squaredExponential','exponential',...
    altitudes,kfgpTimeStep);

kfgp.spatialCovAmp       = 1;
kfgp.spatialLengthScale  = 20;
kfgp.temporalCovAmp      = 1;
kfgp.temporalLengthScale = 10;
kfgp.noiseVariance       = 1e-3;

kfgp.initVals = kfgp.initializeKFGP;
kfgp.spatialCovMat = kfgp.makeSpatialCovarianceMatrix(altitudes);
kfgp.spatialCovMatRoot = kfgp.calcSpatialCovMatRoot;

%% generate synthetic flow data
% number of altitudes
nAlt = numel(altitudes);
% final time for data generation in minutes
tFinData = 600;
% time step for synthetic data generation
timeStepSynData = 1;
% standard deviation for synthetic data generation
stdDevSynData = 0.5;
% get the time series object
[synFlow,synAlt] = kfgp.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);

%% regression using KFGP
% sampling time step
dt = kfgp.kfgpTimeStep;
% algorithm final time
algFinTime = 120;
% sampling time vector
tSamp = 0:dt:algFinTime;
% number of samples
nSamp = numel(tSamp);

% preallocat sampling matrices
xSamp  = NaN(1,nSamp);
ySamp  = NaN(nSamp,1);
flowVal = NaN(nSamp,1);
XTSamp = NaN(2,nSamp);

% preallocate matrices for KFGP
predMeansKFGP = NaN(nAlt,nSamp);
postVarsKFGP  = NaN(nAlt,nSamp);
stdDevKFGP    = NaN(nAlt,nSamp);
upBoundKFGP   = NaN(nAlt,nSamp);
loBoundKFGP   = NaN(nAlt,nSamp);

% number of std deviations for bounds calculations
numStdDev = 1;


%% initialize MPC KFGP
% mpc time step
mpckfgpTimeStep = 1;
% mpc prediction horizon
predictionHorz  = 4;
% fmincon options
options = optimoptions('fmincon','algorithm','sqp','display','off');
% make new KFGP to maintain MPC calculations
mpckfgp = GP.KalmanFilteredGaussianProcess('squaredExponential','exponential',...
    altitudes,mpckfgpTimeStep);

mpckfgp.spatialCovAmp       = kfgp.spatialCovAmp;
mpckfgp.spatialLengthScale  = kfgp.spatialLengthScale;
mpckfgp.temporalCovAmp      = kfgp.temporalCovAmp;
mpckfgp.temporalLengthScale = kfgp.temporalLengthScale;
mpckfgp.noiseVariance       = kfgp.noiseVariance;

mpckfgp.initVals = mpckfgp.initializeKFGP;
mpckfgp.spatialCovMat = mpckfgp.makeSpatialCovarianceMatrix(altitudes);
mpckfgp.spatialCovMatRoot = mpckfgp.calcSpatialCovMatRoot;

mpckfgp.tetherLength         = 100;

% acquistion function parameters
mpckfgp.exploitationConstant = 0;
mpckfgp.explorationConstant  = 2;
mpckfgp.predictionHorizon    = predictionHorz;

% max mean elevation angle step size
uMax = 5;
Astep = zeros(predictionHorz-1,predictionHorz);
bstep = uMax*ones(2*(predictionHorz-1),1);
for ii = 1:predictionHorz-1
    for jj = 1:predictionHorz
        if ii == jj
            Astep(ii,jj) = -1;
            Astep(ii,jj+1) = 1;
        end
        
    end
end
Astep = [Astep;-Astep];
% bounds on first step
fsBoundsA = zeros(2,predictionHorz);
fsBoundsA(1,1) = 1;
fsBoundsA(2,1) = -1;
A = [fsBoundsA;Astep];
% upper and lower bounds
minElev = asin(min(altitudes)/mpckfgp.tetherLength)*180/pi;
lb = minElev*ones(1,predictionHorz);
maxElev = asin(max(altitudes)/mpckfgp.tetherLength)*180/pi;
ub = maxElev*ones(1,predictionHorz);

uAllowable = linspace(-uMax,uMax,5);

jj = 1;

%% do the regresson
for ii = 1:nSamp
    if ii == 1
        nextPoint = altitudes(randperm(nAlt,1));
    end
    
    % go to xSamp
    xSamp(ii) = nextPoint;
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
    else
        % KFGP: initial state estimate
        sk_k   = skp1_kp1;
        ck_k   = ckp1_kp1;
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
    
    % MPCKFGP: make decision about next point
    meanElevation = asin(xSamp(ii)/mpckfgp.tetherLength)*180/pi;
    fsBoundsB(1,1) = meanElevation + uMax;
    fsBoundsB(2,1) = meanElevation - uMax;
    b = [fsBoundsB;bstep];
    
    if mod(tSamp(ii),mpckfgp.kfgpTimeStep)==0
        if ii == 1
            sk_k_mpc = mpckfgp.initVals.s0;
            ck_k_mpc = mpckfgp.initVals.sig0Mat;
        else
            sk_k_mpc = skp1_kp1_mpc;
            ck_k_mpc = ckp1_kp1_mpc;
        end
        
    % KFGP: calculate kalman states
    [F_t_mpc,sigF_t_mpc,skp1_kp1_mpc,ckp1_kp1_mpc] = ...
        mpckfgp.calcKalmanStateEstimates(sk_k,ck_k,xSamp(ii),ySamp(ii));
    % KFGP: calculate prediction mean and posterior variance
    [muKFGP_mpc,sigKFGP_mpc] = mpckfgp.calcPredMeanAndPostVar(altitudes,F_t_mpc,sigF_t_mpc);
    
    F_t_comp = F_t_mpc - F_t;
    sigF_t_comp = sigF_t_mpc - sigF_t;
    sk_k_comp = sk_k_mpc - sk_k;
    ck_k_comp = ck_k_mpc - ck_k;
    skp1_kp1_comp = skp1_kp1_mpc - skp1_kp1;
    ckp1_kp1_comp = ckp1_kp1_mpc - ckp1_kp1;
    muKFGP_comp = muKFGP_mpc - muKFGP;
    sigKFGP_comp = sigKFGP_mpc - sigKFGP;
    
    % KFGP: store them
    predMeansKFGP_MPC(:,jj) = muKFGP_mpc;
    postVarsKFGP_MPC(:,jj)  = sigKFGP_mpc;
    % KFGP: calculate bounds
    stdDevKFGP_MPC(:,jj) = postVarsKFGP_MPC(:,jj).^0.5;
    % KFGP: upper bounds = mean + x*(standard deviation)
    upBoundKFGP_MPC(:,jj) = predMeansKFGP_MPC(:,jj) + numStdDev*stdDevKFGP_MPC(:,jj);
    % KFGP: lower bounds = mean - x*(standard deviation)
    loBoundKFGP_MPC(:,jj) = predMeansKFGP_MPC(:,jj) - numStdDev*stdDevKFGP_MPC(:,jj);

    tMPC(jj) = tSamp(ii);
    jj = jj+1;
    nextPoint = altitudes(randperm(nAlt,1));

    end
    
    
end


%% convert results to time series and store in strcut
regressionRes(1).predMean  = timeseries(predMeansKFGP,tSamp*60);
regressionRes(1).loBound   = timeseries(loBoundKFGP,tSamp*60);
regressionRes(1).upBound   = timeseries(upBoundKFGP,tSamp*60);
regressionRes(1).dataSamp  = timeseries([xSamp;ySamp'],tSamp*60);
regressionRes(1).legend    = 'KFGP';

regressionRes(2).predMean  = timeseries(predMeansKFGP_MPC,tMPC*60);
regressionRes(2).loBound   = timeseries(loBoundKFGP_MPC,tMPC*60);
regressionRes(2).upBound   = timeseries(upBoundKFGP_MPC,tMPC*60);
regressionRes(2).dataSamp  = timeseries([xSamp;ySamp'],tSamp*60);
regressionRes(2).legend    = 'KFGP-MPC';


%% plot the data
F = animatedPlot(synFlow,synAlt,'plotTimeStep',0.05,...
    'regressionResults',regressionRes,'waitforbutton',true);


