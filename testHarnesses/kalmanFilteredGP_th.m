clear
clc
close all

cd(fileparts(mfilename('fullpath')));


%% initialize KFGP
rng(60);

% altitudes
altitudes = 0:10:100;
kfgpTimeStep = 0.05;

% spatial kernel

kfgp = GP.KalmanFilteredGaussianProcess('squaredExponential','exponential',...
    'windPowerLaw',altitudes,kfgpTimeStep);

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
tFinData = 300;
% time step for synthetic data generation
timeStepSynData = 0.5;
% standard deviation for synthetic data generation
stdDevSynData = 0.5;
% get the time series object
[synFlow,synAlt] = kfgp.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);

%% regression using KFGP
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
mpckfgpTimeStep = 0.5;
% mpc prediction horizon
predictionHorz  = 6;
% fmincon options
options = optimoptions('fmincon','algorithm','sqp','display','off');
% make new KFGP to maintain MPC calculations
mpckfgp = GP.KalmanFilteredGaussianProcess('squaredExponential','exponential',...
    'windPowerLaw',altitudes,mpckfgpTimeStep);

mpckfgp.spatialCovAmp       = kfgp.spatialCovAmp;
mpckfgp.spatialLengthScale  = kfgp.spatialLengthScale;
mpckfgp.temporalCovAmp      = kfgp.temporalCovAmp;
mpckfgp.temporalLengthScale = kfgp.temporalLengthScale;
mpckfgp.noiseVariance       = kfgp.noiseVariance;

mpckfgp.initVals = mpckfgp.initializeKFGP;
mpckfgp.spatialCovMat = mpckfgp.makeSpatialCovarianceMatrix(altitudes);
mpckfgp.spatialCovMatRoot = mpckfgp.calcSpatialCovMatRoot;

mpckfgp.tetherLength         = 200;

% acquistion function parameters
mpckfgp.exploitationConstant = 1;
mpckfgp.explorationConstant  = 250;
mpckfgp.predictionHorizon    = predictionHorz;

% max mean elevation angle step size
duMax = 5;
Astep = zeros(predictionHorz-1,predictionHorz);
bstep = duMax*ones(2*(predictionHorz-1),1);
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
minElev = asin(altitudes(2)/mpckfgp.tetherLength)*180/pi;
lb = minElev*ones(1,predictionHorz);
maxElev = asin(max(altitudes)/mpckfgp.tetherLength)*180/pi;
ub = maxElev*ones(1,predictionHorz);

uAllowable = linspace(-duMax,duMax,5);

% number of times mpc will trigger
nMPC = floor(tSamp(end)/mpckfgpTimeStep);
tMPC     = nan(1,nMPC);
jObjFmin = nan(1,nMPC);
jExploitFmin = nan(1,nMPC);
jExploreFmin = nan(1,nMPC);
uTrajFmin = nan(predictionHorz,nMPC);

jObjBF   = nan(1,nMPC);
jExploitBF = nan(1,nMPC);
jExploreBF = nan(1,nMPC);
uTrajBF = nan(predictionHorz,nMPC);

% omniscient controller preallocation
fValOmni     = nan(1,nSamp);
omniElev     = nan(1,nSamp);
elevsAtAllAlts = asin(altitudes/mpckfgp.tetherLength)*180/pi;
cosElevAtAllAlts = cosd(elevsAtAllAlts);
% baseline contoller
baselineElev = 15;
baselineAlt  = mpckfgp.tetherLength*sind(baselineElev);
fValBaseline = nan(1,nSamp);
% KFGP control
fValKFGP     = nan(1,nSamp);
KFGPElev     = nan(1,nSamp);

% mpc counter
jj = 1;
%% do the regresson
for ii = 1:nSamp
    % go to xSamp
    if ii == 1
        nextPoint = altitudes(randperm(nAlt,1));
    end
    xSamp(ii) = nextPoint;
    % measure flow at xSamp(ii) at tSamp(ii)
    fData = resample(synFlow,tSamp(ii)*60).Data;
    hData = resample(synAlt,tSamp(ii)*60).Data;
    ySamp(ii) = interp1(hData,fData,xSamp(ii));
    % calculate pseudo power
    % omniscient, uncontrained controller
    [fValOmni(ii),omniIdx] = max((fData.*cosElevAtAllAlts(:)).^3);
    omniElev(ii) = elevsAtAllAlts(omniIdx);
    % base line
    fValBaseline(ii) = (interp1(hData,fData,baselineAlt)*cosd(baselineElev))^3;
%     % KFGP
    KFGPElev(ii)  = asin(xSamp(ii)/mpckfgp.tetherLength)*180/pi;
    fValKFGP(ii) = (ySamp(ii)*cosd(KFGPElev(ii)))^3;
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
    fsBoundsB(1,1) = meanElevation + duMax;
    fsBoundsB(2,1) = meanElevation - duMax;
    b = [fsBoundsB;bstep];
    
    if ii>1 && mod(tSamp(ii),mpckfgp.kfgpTimeStep)==0
        tMPC(jj) = tSamp(ii);
        % mpc kalman estimate
        [F_t_mpc,sigF_t_mpc,skp1_kp1_mpc,ckp1_kp1_mpc] = ...
            mpckfgp.calcKalmanStateEstimates(sk_k,ck_k,xSamp(ii),ySamp(ii));
        
        % use fminc to solve for best trajectory
        [bestTraj,mpcObj] = ...
            fmincon(@(u) -mpckfgp.calcMpcObjectiveFn(...
            F_t_mpc,sigF_t_mpc,skp1_kp1_mpc,ckp1_kp1_mpc...
            ,u),meanElevation*ones(predictionHorz,1),A,b,[],[]...
            ,lb,ub,[],options);
        uTrajFmin(:,jj) = mpckfgp.calcDelevTraj(meanElevation,bestTraj);
        % get other values
        [jObjFmin(jj),jExptFmin,jExpreFmin] = ...
            mpckfgp.calcMpcObjectiveFn(...
            F_t_mpc,sigF_t_mpc,skp1_kp1_mpc,ckp1_kp1_mpc,...
            bestTraj);
        jExploitFmin(jj) = sum(jExptFmin);
        jExploreFmin(jj) = sum(jExpreFmin);
        
        % brute force
        bruteForceTraj = mpckfgp.bruteForceTrajectoryOpt(...
            F_t_mpc,sigF_t_mpc,skp1_kp1_mpc,ckp1_kp1_mpc,...
            meanElevation,uAllowable,lb(1),ub(1));
        uTrajBF(:,jj) = mpckfgp.calcDelevTraj(meanElevation,bruteForceTraj);
        % get other values
        [jObjBF(jj),jExptBF,jExprBF] = ...
            mpckfgp.calcMpcObjectiveFn(...
            F_t_mpc,sigF_t_mpc,skp1_kp1_mpc,ckp1_kp1_mpc,...
            bruteForceTraj);
        jExploitBF(jj) = sum(jExptBF);
        jExploreBF(jj) = sum(jExprBF);

        % next point
        nextPoint = mpckfgp.convertMeanElevToAlt(bestTraj(1));
%         nextPoint = mpckfgp.convertMeanElevToAlt(bruteForceTraj(1));
        disp(['FMINCON :',num2str(mpckfgp.convertMeanElevToAlt(bestTraj'),'%.3f ')]);
        disp(['BF      :',num2str(mpckfgp.convertMeanElevToAlt(bruteForceTraj),'%.3f ')]);
        fprintf('\n');
        jj = jj+1;
    end
    
    
end


%% convert results to time series and store in strcut
regressionRes(1).predMean  = timeseries(predMeansKFGP,tSamp*60);
regressionRes(1).loBound   = timeseries(loBoundKFGP,tSamp*60);
regressionRes(1).upBound   = timeseries(upBoundKFGP,tSamp*60);
regressionRes(1).dataSamp  = timeseries([xSamp;ySamp'],tSamp*60);
regressionRes(1).dataAlts  = synAlt;
regressionRes(1).legend    = 'KFGP';

meanVals = [mean(fValOmni) mean(fValBaseline) mean(fValKFGP)];

fName = ['results\KFGPres ',strrep(datestr(datetime),':','-'),'.mat'];
save(fName);

%% plot the data
% look at objective function values at each step of MPC
figure

spOpt = gobjects;
spOpt(1) = subplot(1,3,1);
plot(1:nMPC,jObjBF,'-o')
grid on;hold on
plot(1:nMPC,jObjFmin,'-o')
legend('BF','FMINCON')
xlabel('MPC step');
ylabel('$J_{total}$')


spOpt(2) = subplot(1,3,2);
plot(1:nMPC,jExploitBF,'-o')
grid on;hold on
plot(1:nMPC,jExploitFmin,'-o')
legend('BF','FMINCON')
xlabel('MPC step');
ylabel('$J_{exploit}$')


spOpt(3) = subplot(1,3,3);
plot(1:nMPC,jExploreBF,'-o')
grid on;hold on
plot(1:nMPC,jExploreFmin,'-o')
legend('BF','FMINCON')
xlabel('MPC step');
ylabel('$J_{explore}$')

linkaxes(spOpt(1:3),'x')

sgtitle('Optimization algorithms comparisons')


% fval plots
fvalPlot = gobjects;
legendStrs = {'Omniscient unconstrained','Baseline','KFGP constrained'};

figure
subplot(1,2,1)
fvalPlot(1) = plot(0:nSamp-1,fValOmni,'linewidth',1);
grid on;hold on
fvalPlot(2) = plot(0:nSamp-1,fValBaseline,'linewidth',1);
fvalPlot(3) = plot(0:nSamp-1,fValKFGP,'linewidth',1);
legend(fvalPlot(1:3),legendStrs);
xlabel('Time step');
ylabel('$(v_{f} cos(\phi))^3$')


subplot(1,2,2)
elevPlots = gobjects;
elevPlots(1) = stairs(0:nSamp-1,omniElev,'linewidth',1);
grid on;hold on
elevPlots(2) = plot([0 nSamp-1],baselineElev*[1 1],'linewidth',1);
elevPlots(3) = stairs(0:nSamp-1,KFGPElev,'linewidth',1);
legend(elevPlots(1:3),legendStrs);
xlabel('Time step');
ylabel('$\phi$')
sgtitle(sprintf('MPC triggered every %d steps',mpckfgpTimeStep/kfgpTimeStep));

set(findobj('-property','FontSize'),'FontSize',11);

%% animation
figure
F = animatedPlot(synFlow,synAlt,'plotTimeStep',0.25,...
    'regressionResults',regressionRes...
    ,'waitforbutton',true);


