clear
clc
% close all

cd(fileparts(mfilename('fullpath')));
% load('envFile');


%% initialize GP
rng(5);

gp = GP.GaussianProcess('squaredExponential','squaredExponential','windPowerLaw');

gp.spatialCovAmp       = 5.1^2;
gp.spatialLengthScale  = 220;
gp.temporalCovAmp      = 1;
gp.temporalLengthScale = 22;
gp.noiseVariance       = 1e-3;

%% generate synthetic flow data
% altitudes
altitudes = 0:50:1000;
nAlt = numel(altitudes);
% final time for data generation in minutes
tFinData = 300;
% time step for synthetic data generation
timeStepSynData = 1;
% standard deviation for synthetic data generation
stdDevSynData = 4;
% get the time series object
[synFlow,synAlt] = gp.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);

%% process data
figure;
axisObj = paperFlowPlot(synFlow,altitudes,30,180);

%% export file
saveFile = input('Save file? Options: Enter y or n\n','s');
if strcmpi(saveFile,'y')
filName = strcat('flowImage_',strrep(datestr(datetime),':','-'));
savefig(filName);
exportgraphics(axisObj,[filName,'.png'],'Resolution',600)
end
