clear
clc
close all

cd(fileparts(mfilename('fullpath')));
load('envFile');


%% initialize GP
rng(1);

gp = GP.GaussianProcess('squaredExponential','exponential','windPowerLaw');

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
stdDevSynData = 2;
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
% F = animatedPlot(synFlow,synAlt,'plotTimeStep',1);

%% test CNAPS flow
cnapsData.dataHours = 10;
cnapsData.dataIdx   = testEnv.flowVecTimeseries.Value.Time<=cnapsData.dataHours*3600;
cnapsTvec    = squeeze(testEnv.flowVecTimeseries.Value.Time(cnapsData.dataIdx));
cnapsFlow    = squeeze(testEnv.flowVecTimeseries.Value.Data(1,1,:,1,cnapsData.dataIdx));
synCnaps     = timeseries(cnapsFlow,cnapsTvec);

CNAPSalts    = timeseries(repmat(testEnv.zGridPoints.Value',1,2),...
    [cnapsTvec(1),cnapsTvec(end)]);

F2 = animatedPlot(synCnaps,CNAPSalts,'plotTimeStep',10,'wait',false);

%% video settings

