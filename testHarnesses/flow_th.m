clear
clc
close all

cd(fileparts(mfilename('fullpath')));


%% initialize GP
rng(1);

gp = GP.GaussianProcess('squaredExponential','exponential');

gp.spatialCovAmp       = 1;
gp.spatialLengthScale  = 30;
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


%% generate flow data as a function of elevation angle
% tether length
thrLength = 200;
% elevation angles
elevAngles = asin((altitudes./thrLength).^1)*180/pi;

gp2 = gp;
gp2.spatialLengthScale = asin((gp.spatialLengthScale/thrLength)^0.5)*180/pi;

% get the time series object
[synFlow2,synElev] = gp2.generateSyntheticFlowData(elevAngles,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);

%% plot the data
F = animatedPlot(synFlow,synAlt,'plotTimeStep',1);

