clear
clc
format compact
close all

%% change working directory to current folder
cd(fileparts(mfilename('fullpath')));

%% generate wind using colored noise
% rng(56);
% environment
hMax = 1500;
hMin = 100;
heights = hMin:100:hMax;
heights = heights(:);
meanFlow = 0;
noTP = numel(heights);
% time in minutes
timeStep = 0.05*5;
tVec = 0:timeStep:1*60;
noTimeSteps = numel(tVec);
% time in seconds
timeInSec = 60*tVec;
% std deviation for wind data generation
stdDev = 1;
% hyper parameters
timeScale = 10;
heightScale = 200;
% generate data
windSpeedOut = meanFlow + genWindv2(heights,heightScale,tVec,timeScale,stdDev);
heights2 = repmat(heights,1,noTimeSteps);
tVec2 = repmat(tVec(:)',noTP,1);
dsgnPts = [heights2(:)'; tVec2(:)'];
dsgnFvals = windSpeedOut(:);

% translate and scale wind data
Flows = NaN(noTP,1,noTimeSteps);
for ii = 1:noTimeSteps
    Flows(:,:,ii) = windSpeedOut(:,ii);
end

simTime = 1*60*60;
Heights = timeseries(repmat(heights,1,1,2),[0 simTime]);
Flows = timeseries(Flows,timeInSec);

%% set up the KFGP classdef
% construct an instance of the RGP class: squaredExponential or exponential
% acquisition function options: expectedImprovement or upperConfidenceBound
EstimationGpkf = GPKF(1,'exponential','upperConfidenceBound');
% set values of hyper parameters
EstimationGpkf.p_noiseVariance = 0.01;
EstimationGpkf.p_spatialCovarianceAmp = 1;
EstimationGpkf.p_spatialLenghtScale = heightScale;
EstimationGpkf.p_temporalLenghtScale = timeScale;

% set gpkf parameters
EstimationGpkf.p_gpfkTimeStep = timeStep;
EstimationGpkf.p_xMeasure = heights';
EstimationGpkf.p_squaredExpApproxOrder = 2;

% set acquisition function parameters
EstimationGpkf.p_exploitationConstant = 0;
EstimationGpkf.p_explorationConstant = 4;

% set initial values
EstimationGpkf = EstimationGpkf.m_CalcSpatialCovMat;
EstimationGpkf = EstimationGpkf.m_CalcSpatialCovRoot;
EstimationGpkf = EstimationGpkf.gpfkInitialize;

%% set up MPC GPKF
% % % initialize MPC GPKF to be the same as estimation gpkf
MpcGpkf = EstimationGpkf;
% % % change the sampling time step
MpcGpkf.p_gpfkTimeStep = timeScale/10;
% set initial values
MpcGpkf = MpcGpkf.m_CalcSpatialCovMat;
MpcGpkf = MpcGpkf.m_CalcSpatialCovRoot;
MpcGpkf = MpcGpkf.gpfkInitialize;
% % % prediction horizon
predHorz = 5;
% % % allowable control inputs
uStep = 100;
uAllowable = uStep*(-2:1:2);

%% do the actual gaussian process kalman filtering
% % % make a finer domain over which predictions are made
xPredict = linspace(heights(1),heights(end),1*numel(heights));

ck_k = EstimationGpkf.p_gpfkInitVals.sig0Mat;
sk_k = EstimationGpkf.p_gpfkInitVals.s0;

% % % number of iterations
noIter = noTimeSteps;
% % % preallocate matrices
predMean = NaN(size(xPredict,2),noIter);
postVar =  predMean;
acquiFun = predMean;
stdDev =  predMean;
upperBound = predMean;
lowerBound = predMean;
pointsVisited = NaN(1,noIter);
fValAtPt = NaN(noIter,1);
% start off mpc counter
mpcCount = 1;
numStarts = 4;

for ii = 1:noIter
    % % % visit said points
    if ii == 1 || mpcCount == 1
        visitIdx = randperm(size(EstimationGpkf.p_xMeasure,2),1);
        Mk = EstimationGpkf.p_xMeasure(:,visitIdx);
    else
        [~,maxIdx] =  max(postVar(:,ii-1));
        Mk2 = EstimationGpkf.p_xMeasure(:,maxIdx);
        [~,visitIdx] = min(abs(chosenTrajectory(1)-...
            EstimationGpkf.p_xMeasure));
        Mk = EstimationGpkf.p_xMeasure(visitIdx);
        
    end
    % % % extract wind speed at visited values
    yk = windSpeedOut((Mk == EstimationGpkf.p_xMeasure)',ii);
    % % % gpkf MPC
    if mod(tVec(ii),MpcGpkf.p_gpfkTimeStep) == 0 && tVec(ii)~=0
        
        
        mpcFmincon = MpcGpkf.m_gpkfMPC_fmincon(sk_k,ck_k,Mk,yk,...
            uAllowable,predHorz,numStarts);
        
        mpcResBF = MpcGpkf.m_gpkfMPC_bruteForce(sk_k,ck_k,Mk,yk,...
            uAllowable,predHorz);

        
        % % % store optimal control sequence)
        optStateTrajBF(mpcCount,:) = mpcResBF.optStateTrajectory;
        optFvalBF(mpcCount) = mpcResBF.objFunVal;
        optStateTrajFmin(mpcCount,:) = mpcFmincon.optStateTrajectoryFmin;
        optFvalFmin(mpcCount) = mpcFmincon.objFunValFmin;
        
        chosenTrajectory = mpcFmincon.optStateTrajectoryFmin;
%         chosenTrajectory = mpcResBF.optStateTrajectory;
        mpcCount = mpcCount + 1;
    end
    % % % kalman state estimation
    [F_t,sigF_t,skp1_kp1,ckp1_kp1] = EstimationGpkf.m_gpkfKalmanEstimation...
        (sk_k,ck_k,Mk,yk);
    
    % % % regression over a finer domain
    [predMean(:,ii),postVar(:,ii)] = EstimationGpkf.m_gpkfRegression(...
        xPredict,F_t,sigF_t);
    
    % % % remove real or imaginary parts lower than eps
    stdDev(:,ii) = postVar(:,ii).^0.5;
    % % % upper bounds = mean + x*(standard deviation)
    upperBound(:,ii) = predMean(:,ii) + 1*stdDev(:,ii);
    % % % lower bounds = mean + x*(standard deviation)
    lowerBound(:,ii) = predMean(:,ii) - 1*stdDev(:,ii);
    % % % calculate acquisition function value
    acquiFun(:,ii) = EstimationGpkf.m_calcAcquisitionFun(...
        predMean(:,ii),postVar(:,ii),yk);
    
    % % % update previous step information
    sk_k = skp1_kp1;
    ck_k = ckp1_kp1;
    
    % % % store points visited at the respective function value
    pointsVisited(:,ii) = Mk;
    fValAtPt(ii,:) = yk;
    
end

%% save data to output folder
[status, msg, msgID] = mkdir(pwd,'outputs');
fName = [pwd,'\outputs\',strrep(datestr(datetime),':','_')];

% delete([pwd,'\outputs\*.mat'])
% delete([pwd,'\outputs\*.avi'])

save(fName)

%% plot file
makePlotsAndAnimations


