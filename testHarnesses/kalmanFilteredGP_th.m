clear
clc
close all

cd(fileparts(mfilename('fullpath')));


%% initialize GP
rng(1);

kfgp = GP.KalmanFilteredGaussianProcess(1,'squaredExponential','exponential');

kfgp.spatialCovAmp       = 1;
kfgp.spatialLengthScale  = 20;
kfgp.temporalCovAmp      = 1;
kfgp.temporalLengthScale = 10;
kfgp.noiseVariance       = 1e-3;

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
[synFlow,synAlt] = kfgp.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);


