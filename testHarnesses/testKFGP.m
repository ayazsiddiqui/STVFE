clear;clc;
% close all;
cd(fileparts(mfilename('fullpath')));

fIdx = 1;

%% initailize
% load vehicle
load('ampyxVhcl.mat');

% initialize class
cIn = maneuverabilityAdvanced(vhcl);
cIn.wingOswaldEff   = 0.3;
cIn.hstabOswaldEff  = 0.6;
cIn.hstabZerAoADrag = 0.1*cIn.hstabZerAoADrag;
cIn.vstabOswaldEff  = 0.3;

% tether length
cIn.tetherLength = 1000;
cIn.pathWidth    = 20;
cIn.pathHeight   = 12;

% turbine properties
cIn.turbCP = 0.5;
cIn.turbCD = 1.5*cIn.turbCP;
% calculate optimum turbine diameter
cIn.turbDia = cIn.calcOptTurbDiameter(11*pi/180);

% path parameter and some other constants
S.nPoints = 41;
S.pathParam = linspace(0,2*pi,S.nPoints);

% tangent pitch angle
S.tgtPitch = 6*pi/180;
S.elevatorDeflection = 0;

%% initialize KFGP
rng(8);

% altitudes
altitudes = 0:100:1000;
kfgpTimeStep = 0.2;

% spatial kernel
spaceKernel = 'squaredExponential';
timeKernel  = 'squaredExponential';

kfgp = GP.KalmanFilteredGaussianProcess(spaceKernel,timeKernel,...
    'windPowerLaw',altitudes,kfgpTimeStep);

kfgp.spatialCovAmp       = 5.1^2;
kfgp.spatialLengthScale  = 220;
kfgp.temporalCovAmp      = 1;
kfgp.temporalLengthScale = 22;
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
timeStepSynData = 3;
% standard deviation for synthetic data generation
stdDevSynData = 4;
% get the time series object
[synFlow,synAlt] = kfgp.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);

%% regression using KFGP
% algorithm final time
algFinTime = 180;
% sampling time vector
tSamp = 0:kfgp.kfgpTimeStep:algFinTime;
% number of samples
nSamp = numel(tSamp);

% preallocat sampling matrices
xSamp   = NaN(1,nSamp);
ySamp   = NaN(nSamp,1);
flowVal = NaN(nSamp,1);
XTSamp  = NaN(2,nSamp);

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
mpckfgpTimeStep = 3;
% mpc prediction horizon
predictionHorz  = 6;
% fmincon options
options = optimoptions('fmincon','algorithm','sqp','display','off');
% make new KFGP to maintain MPC calculations
mpckfgp = GP.KalmanFilteredGaussianProcess(spaceKernel,timeKernel,...
    'windPowerLaw',altitudes,mpckfgpTimeStep);

mpckfgp.spatialCovAmp       = kfgp.spatialCovAmp;
mpckfgp.spatialLengthScale  = kfgp.spatialLengthScale;
mpckfgp.temporalCovAmp      = kfgp.temporalCovAmp;
mpckfgp.temporalLengthScale = kfgp.temporalLengthScale;
mpckfgp.noiseVariance       = kfgp.noiseVariance;

mpckfgp.initVals            = mpckfgp.initializeKFGP;
mpckfgp.spatialCovMat       = mpckfgp.makeSpatialCovarianceMatrix(altitudes);
mpckfgp.spatialCovMatRoot   = mpckfgp.calcSpatialCovMatRoot;

mpckfgp.tetherLength        = cIn.tetherLength;

% acquistion function parameters
mpckfgp.exploitationConstant = 1;
mpckfgp.explorationConstant  = 1;
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
minElev = 5;
lb      = minElev*ones(1,predictionHorz);
maxElev = 60;
ub      = maxElev*ones(1,predictionHorz);

uAllowable = linspace(-duMax,duMax,5);

% number of times mpc will trigger
nMPC         = floor(tSamp(end)/mpckfgpTimeStep);
tMPC         = nan(1,nMPC);
jObjFmin     = nan(1,nMPC);
jExploitFmin = nan(1,nMPC);
jExploreFmin = nan(1,nMPC);
uTrajFmin    = nan(predictionHorz,nMPC);

jObjBF     = nan(1,nMPC);
jExploitBF = nan(1,nMPC);
jExploreBF = nan(1,nMPC);
uTrajBF    = nan(predictionHorz,nMPC);

% omniscient controller preallocation
fValOmni       = nan(1,nSamp);
runAvgOmni     = nan(1,nSamp);
omniElev       = nan(1,nSamp);
elevsAtAllAlts = min(max(minElev,asin(altitudes/mpckfgp.tetherLength)*180/pi),maxElev);
omniAlts = mpckfgp.convertMeanElevToAlt(elevsAtAllAlts);
cosElevAtAllAlts = cosd(elevsAtAllAlts);
meanFnVec = kfgp.meanFunction(altitudes);

% baseline contoller
fValBaseline   = nan(1,nSamp);
runAvgbaseline = nan(1,nSamp);
% KFGP control
fValKFGP     = nan(1,nSamp);
KFGPElev     = nan(1,nSamp);
runAvgKFGP   = nan(1,nSamp);

% path speed analysis results
SR.pathSpeeds     = nan(S.nPoints,nSamp);
SR.pathRoll       = nan(S.nPoints,nSamp);
SR.avgSpeed       = nan(1,nSamp);
SR.avgV_appxCubed = nan(1,nSamp);
SR.avgPower       = nan(1,nSamp);
SR.lapTime        = nan(1,nSamp);
SR.avgAoA         = nan(1,nSamp);
SR.maxTangentRoll = nan(1,nSamp);
% mpc counter
jj = 1;

%% omniscient
SRO = SR;
for ii = 1:nSamp
    % measure flow at xSamp(ii) at tSamp(ii)
    fData = resample(synFlow,tSamp(ii)*60).Data;
    hData = resample(synAlt,tSamp(ii)*60).Data;
    % calculate pseudo power
    % omniscient, uncontrained controller
    omnifData = interp1(hData,fData,omniAlts);
    [fValOmni(ii),omniIdx] = max(cosineFlowCubed(omnifData,cosElevAtAllAlts));
    runAvgOmni(ii) = mean(fValOmni(1:ii));
    omniElev(ii) = elevsAtAllAlts(omniIdx);
    
    cIn.meanElevationInRadians = omniElev(ii)*pi/180;
    solVals = cIn.getAttainableVelocityOverPath(interp1(hData,fData,omniAlts(omniIdx)),...
        S.tgtPitch,S.pathParam);
    % calc relevant results
    SRO.pathSpeeds(:,ii)   = solVals.vH_path;
    SRO.pathRoll(:,ii)     = solVals.roll_path;
    SRO.avgSpeed(ii)       = mean(solVals.vH_path);
    SRO.avgV_appxCubed(ii) = mean(max(0,-solVals.B_Vapp_path(1,:)).^3);
    SRO.avgPower(ii)       = cIn.turbCP*0.5*cIn.fluidDensity*...
        SRO.avgV_appxCubed(ii)*cIn.turbArea;
    SRO.lapTime(ii)        = cIn.pathLength/SRO.avgSpeed(ii);
    SRO.avgAoA(ii)         = mean(cIn.calcAngleOfAttackInRadians(...
        solVals.B_Vapp_path))*180/pi;
    SRO.maxTangentRoll(ii) = max(abs(SRO.pathRoll(:,ii)))*180/pi;
end

%% baseline
baselineElev   = ceil(mean(omniElev));
baselineAlt    = mpckfgp.tetherLength*sind(baselineElev);

SRB = SR;
for ii = 1:nSamp
    % measure flow at xSamp(ii) at tSamp(ii)
    fData = resample(synFlow,tSamp(ii)*60).Data;
    hData = resample(synAlt,tSamp(ii)*60).Data;
    % base line
    fValBaseline(ii) = cosineFlowCubed(interp1(hData,fData,baselineAlt),cosd(baselineElev));
    runAvgbaseline(ii) = mean(fValBaseline(1:ii));
    
    cIn.meanElevationInRadians = baselineElev*pi/180;
    solVals = cIn.getAttainableVelocityOverPath(interp1(hData,fData,baselineAlt),...
        S.tgtPitch,S.pathParam);
    % calc relevant results
    SRB.pathSpeeds(:,ii)   = solVals.vH_path;
    SRB.pathRoll(:,ii)     = solVals.roll_path;
    SRB.avgSpeed(ii)       = mean(solVals.vH_path);
    SRB.avgV_appxCubed(ii) = mean(max(0,-solVals.B_Vapp_path(1,:)).^3);
    SRB.avgPower(ii)       = cIn.turbCP*0.5*cIn.fluidDensity*...
        SRB.avgV_appxCubed(ii)*cIn.turbArea;
    SRB.lapTime(ii)        = cIn.pathLength/SRB.avgSpeed(ii);
    SRB.avgAoA(ii)         = mean(cIn.calcAngleOfAttackInRadians(...
        solVals.B_Vapp_path))*180/pi;
    SRB.maxTangentRoll(ii) = max(abs(SRB.pathRoll(:,ii)))*180/pi;
end

%% do the regresson
for ii = 1:nSamp
    % go to xSamp
    if ii == 1
        nextPoint = baselineAlt;
    end
    xSamp(ii) = nextPoint;
    % measure flow at xSamp(ii) at tSamp(ii)
    fData = resample(synFlow,tSamp(ii)*60).Data;
    hData = resample(synAlt,tSamp(ii)*60).Data;
    flowVal(ii) = interp1(hData,fData,xSamp(ii));
    ySamp(ii) =  kfgp.meanFunction(xSamp(ii)) - flowVal(ii);
    % calculate pseudo power
    % KFGP
    KFGPElev(ii)  = asin(xSamp(ii)/mpckfgp.tetherLength)*180/pi;
    fValKFGP(ii)  = cosineFlowCubed(flowVal(ii),cosd(KFGPElev(ii)));
    runAvgKFGP(ii) = mean(fValKFGP(1:ii));
    cIn.meanElevationInRadians = KFGPElev(ii)*pi/180;
    solVals = cIn.getAttainableVelocityOverPath(flowVal(ii),...
        S.tgtPitch,S.pathParam);
    % calc relevant results
    SR.pathSpeeds(:,ii)   = solVals.vH_path;
    SR.pathRoll(:,ii)     = solVals.roll_path;
    SR.avgSpeed(ii)       = mean(solVals.vH_path);
    SR.avgV_appxCubed(ii) = mean(max(0,-solVals.B_Vapp_path(1,:)).^3);
    SR.avgPower(ii)       = cIn.turbCP*0.5*cIn.fluidDensity*...
        SR.avgV_appxCubed(ii)*cIn.turbArea;
    SR.lapTime(ii)        = cIn.pathLength/SR.avgSpeed(ii);
    SR.avgAoA(ii)         = mean(cIn.calcAngleOfAttackInRadians(...
        solVals.B_Vapp_path))*180/pi;
    SR.maxTangentRoll(ii) = max(abs(SR.pathRoll(:,ii)))*180/pi;
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
    predMeansKFGP(:,ii) = meanFnVec(:) - muKFGP;
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
    fsBoundsB(2,1) = -(meanElevation - duMax);
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
        
        % next point
        nextPoint = mpckfgp.convertMeanElevToAlt(bestTraj(1));
%         nextPoint = mpckfgp.convertMeanElevToAlt(bruteForceTraj(1));
        disp(['FMINCON :',num2str(mpckfgp.convertMeanElevToAlt(bestTraj'),'%.3f ')]);
%         disp(['BF      :',num2str(mpckfgp.convertMeanElevToAlt(bruteForceTraj),'%.3f ')]);
        fprintf('\n');
        jj = jj+1;
    end
    
end

SR.powerRunningAvg = nan(1,nSamp);
SRO.powerRunningAvg = nan(1,nSamp);
SRB.powerRunningAvg = nan(1,nSamp);


for ii = 1:nSamp
    SR.powerRunningAvg(ii) = mean(SR.avgPower(1:ii))/1000;
    SRO.powerRunningAvg(ii) =  mean(SRO.avgPower(1:ii))/1000;
    SRB.powerRunningAvg(ii) =  mean(SRB.avgPower(1:ii))/1000;
    
end


%% convert results to time series and store in strcut
regressionRes(1).predMean  = timeseries(predMeansKFGP,tSamp*60);
regressionRes(1).loBound   = timeseries(loBoundKFGP,tSamp*60);
regressionRes(1).upBound   = timeseries(upBoundKFGP,tSamp*60);
regressionRes(1).dataSamp  = timeseries([xSamp;flowVal'],tSamp*60);
regressionRes(1).dataAlts  = synAlt;
regressionRes(1).legend    = 'KFGP';

meanVals = [mean(fValOmni) mean(fValBaseline) mean(fValKFGP)];

basePerc = 100*mean(fValBaseline)/mean(fValOmni);
mpcPerc = 100*mean(fValKFGP)/mean(fValOmni);
fprintf('du exp baseline mpc\n');
fprintf('%0.0f & %.0f & %.2f & %.2f\n',...
    [duMax mpckfgp.explorationConstant basePerc mpcPerc]);

% fName = ['results\KFGPres ',strrep(datestr(datetime),':','-'),'.mat'];
% save(fName);

%% plot the data
cols = [228,26,28
    77,175,74
    55,126,184]/255;

spIdx = 1;
spAxes = gobjects;
spObj = gobjects;

lwd = 1.2;
tSampPlot = tSamp/60;

% plot elevation angle trajectory
pIdx = 1;
spAxes(spIdx) = subplot(3,1,spIdx);
hold(spAxes(spIdx),'on');
ylabel(spAxes(spIdx),'$\mathbf{\theta_{sp}}$ \textbf{[deg]}','fontweight','bold');
spObj(pIdx) = stairs(tSampPlot,omniElev...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot([tSampPlot(1) tSampPlot(end)],baselineElev*[1 1]...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = stairs(tSampPlot,KFGPElev...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);

% plot instantaenuos jExploit
spIdx = spIdx + 1;
spAxes(spIdx) = subplot(3,1,spIdx);
pIdx = 1;
hold(spAxes(spIdx),'on');
ylabel(spAxes(spIdx),'$\mathbf{J_{exploit}(t_{k})}$','fontweight','bold');
spObj(pIdx) = plot(tSampPlot,fValOmni...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot(tSampPlot,fValBaseline...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot(tSampPlot,fValKFGP...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);

% plot running average j_exploit
spIdx = spIdx + 1;
spAxes(spIdx) = subplot(3,1,spIdx);
pIdx = 1;
hold(spAxes(spIdx),'on');
ylabel(spAxes(spIdx),'\textbf{Avg.} $\mathbf{J_{exploit}}$','fontweight','bold');
spObj(pIdx) = plot(tSampPlot,runAvgOmni...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot(tSampPlot,runAvgbaseline...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot(tSampPlot,runAvgKFGP...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
legend('Omniscient','Baseline','MPC','location','bestoutside'...
    ,'orientation','horizontal')

% axes props
grid(spAxes(1:end),'on');
set(spAxes(1:end),'GridLineStyle',':')
xlabel(spAxes(1:end),'\textbf{Time [hr]}','fontweight','bold');   
set(spAxes(1:end),'FontSize',12);
% spAxes(1).YTick = linspace(spAxes(1).YTick(1),spAxes(1).YTick(end),3);
% spAxes(2).YTick = linspace(spAxes(2).YTick(1),spAxes(2).YTick(end),3);
set(gcf,'InnerPosition',1*[-00 -00 560 1.8*420])

spAxes(1).YLabel.Position(1) = spAxes(2).YLabel.Position(1);
spAxes(3).YLabel.Position(1) = spAxes(2).YLabel.Position(1);


figure
plot(tSampPlot,SR.powerRunningAvg)
hold on
plot(tSampPlot,SRO.powerRunningAvg)
plot(tSampPlot,SRB.powerRunningAvg)

figure
plot(tSampPlot,runAvgKFGP/max(runAvgKFGP))
hold on
plot(tSampPlot,SR.powerRunningAvg/max(SR.powerRunningAvg),'--')
grid on
legend('Surrogate','steady flight tool')

%% animation
figure
F = animatedPlot(synFlow,synAlt,'plotTimeStep',kfgpTimeStep,...
    'regressionResults',regressionRes...
    ,'waitforbutton',false);


