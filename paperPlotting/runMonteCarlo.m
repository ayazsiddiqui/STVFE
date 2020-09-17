clear
clc
% close all

cd(fileparts(mfilename('fullpath')));


%% initialize KFGP

% altitudes
altitudes = 0:100:1000;
kfgpTimeStep = 0.05;

% spatial kernel

kfgp = GP.KalmanFilteredGaussianProcess('squaredExponential','exponential',...
    'windPowerLaw',altitudes,kfgpTimeStep);

kfgp.spatialCovAmp       = 5.1^2;
kfgp.spatialLengthScale  = 220;
kfgp.temporalCovAmp      = 1;
kfgp.temporalLengthScale = 20;
kfgp.noiseVariance       = 1e-3;

kfgp.initVals = kfgp.initializeKFGP;
kfgp.spatialCovMat = kfgp.makeSpatialCovarianceMatrix(altitudes);
kfgp.spatialCovMatRoot = kfgp.calcSpatialCovMatRoot;

%% generate synthetic flow data
% number of altitudes
nAlt = numel(altitudes);
% final time for data generation in minutes
tFinData = 120;
% time step for synthetic data generation
timeStepSynData = 1;
% standard deviation for synthetic data generation
stdDevSynData = 4;

%% regression using KFGP
% algorithm final time
algFinTime = 60;
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
mpckfgpTimeStep = 1;
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

mpckfgp.initVals            = mpckfgp.initializeKFGP;
mpckfgp.spatialCovMat       = mpckfgp.makeSpatialCovarianceMatrix(altitudes);
mpckfgp.spatialCovMatRoot   = mpckfgp.calcSpatialCovMatRoot;

mpckfgp.tetherLength        = 1000;

% acquistion function parameters
mpckfgp.exploitationConstant = 1;
mpckfgp.predictionHorizon    = predictionHorz;

% max mean elevation angle step size
duMax = 2;
Astep = zeros(predictionHorz-1,predictionHorz);
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


%% monte carlo setup
nDataSets = 5;
rngSeeds = randi(100,[nDataSets,1]);

% mobility study
duMaxSweep = 2:10:20;
betaSweep  = 0:100:300;

[DUMAX,BETA] = meshgrid(duMaxSweep,betaSweep);
KFGPFVAL     = nan*DUMAX;

moteCarloRes.omniscient = nan(1,1,nDataSets);
moteCarloRes.baseline   = nan(1,1,nDataSets);
moteCarloRes.kfgpMPC    = nan(size(DUMAX,1),size(DUMAX,2),nDataSets);


for cc = 1:nDataSets
    
    rng(rngSeeds(cc));
    
    % get the time series object
    [synFlow,synAlt] = kfgp.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
        'timeStep',timeStepSynData);
    
    
    %% omniscient
    for oo = 1:nSamp
        % measure flow at xSamp(ii) at tSamp(ii)
        fData = resample(synFlow,tSamp(oo)*60).Data;
        hData = resample(synAlt,tSamp(oo)*60).Data;
        % calculate pseudo power
        % omniscient, uncontrained controller
        omnifData = interp1(hData,fData,omniAlts);
        [fValOmni(oo),omniIdx] = max(cosineFlowCubed(omnifData,cosElevAtAllAlts));
        runAvgOmni(oo) = mean(fValOmni(1:oo));
        omniElev(oo) = elevsAtAllAlts(omniIdx);
    end
    
    moteCarloRes.omniscient(cc) = mean(fValOmni);
    
    %% baseline
    baselineElev   = ceil(mean(omniElev));
    baselineAlt    = mpckfgp.tetherLength*sind(baselineElev);
    
    for bb = 1:nSamp
        % measure flow at xSamp(ii) at tSamp(ii)
        fData = resample(synFlow,tSamp(bb)*60).Data;
        hData = resample(synAlt,tSamp(bb)*60).Data;
        % base line
        fValBaseline(bb) = cosineFlowCubed(interp1(hData,fData,baselineAlt),cosd(baselineElev));
        runAvgbaseline(bb) = mean(fValBaseline(1:bb));
    end
    
    moteCarloRes.baseline(cc) = mean(fValBaseline);

    
    %% do the regresson
    for mm = 1:numel(DUMAX)
        bstep = DUMAX(mm)*ones(2*(predictionHorz-1),1);
        mpckfgp.explorationConstant  = BETA(mm);
        
        % mpc counter
        jj = 1;
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
                jj = jj+1;
            end
            
        end
        
        KFGPFVAL(mm) = mean(fValKFGP); 
        KFGP_OMNI    = mean(fValKFGP)/moteCarloRes.omniscient(cc);
        
    end
    
    moteCarloRes.kfgpMPC(:,:,cc) = KFGPFVAL;

end

filName = strcat('monteCarloStruct_',strrep(datestr(datetime),':','-'));
save(filName,'moteCarloRes','rngSeeds','DUMAX','BETA');

%% convert results to time series and store in strcut
AVG_VALS.KFGP = mean(moteCarloRes.kfgpMPC,3);
AVG_VALS.OMNI = mean(moteCarloRes.omniscient);
AVG_VALS.BASE = mean(moteCarloRes.baseline);

PERC_VALS.KFGP_OMNI = mean(moteCarloRes.kfgpMPC./moteCarloRes.omniscient,3);
PERC_VALS.KFGP_BASE = mean(moteCarloRes.kfgpMPC./moteCarloRes.baseline,3);

contourf(DUMAX,BETA,PERC_VALS.KFGP_OMNI);

%% plot the data
cols = [228,26,28
    77,175,74
    55,126,184]/255;

spIdx = 1;
spAxes = gobjects;
spObj = gobjects;

lwd = 1.2;

% plot elevation angle trajectory
pIdx = 1;
spAxes(spIdx) = subplot(3,1,spIdx);
hold(spAxes(spIdx),'on');
ylabel(spAxes(spIdx),'$\mathbf{\gamma_{sp}}$ \textbf{[deg]}','fontweight','bold');
spObj(pIdx) = stairs(tSamp,omniElev...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot([tSamp(1) tSamp(end)],baselineElev*[1 1]...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = stairs(tSamp,KFGPElev...
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
spObj(pIdx) = plot(tSamp,fValOmni...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot(tSamp,fValBaseline...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot(tSamp,fValKFGP...
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
spObj(pIdx) = plot(tSamp,runAvgOmni...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot(tSamp,runAvgbaseline...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
pIdx = pIdx + 1;
spObj(pIdx) = plot(tSamp,runAvgKFGP...
    ,'-'...
    ,'color',cols(pIdx,:)...
    ,'MarkerFaceColor',cols(pIdx,:)...
    ,'linewidth',lwd);
legend('Omniscient','Baseline','MPC','location','bestoutside'...
    ,'orientation','horizontal')

% axes props
grid(spAxes(1:end),'on');
set(spAxes(1:end),'GridLineStyle',':')
xlabel(spAxes(1:end),'\textbf{Time [min]}','fontweight','bold');
set(spAxes(1:end),'FontSize',12);
% spAxes(1).YTick = linspace(spAxes(1).YTick(1),spAxes(1).YTick(end),3);
% spAxes(2).YTick = linspace(spAxes(2).YTick(1),spAxes(2).YTick(end),3);
set(gcf,'InnerPosition',1*[-00 -00 560 1.8*420])

spAxes(1).YLabel.Position(1) = spAxes(2).YLabel.Position(1);
spAxes(3).YLabel.Position(1) = spAxes(2).YLabel.Position(1);

%%
% saveFile = input('Save file? Options: Enter y or n\n','s');
% if strcmpi(saveFile,'y')
% filName = strcat('ctrlRes_',strrep(datestr(datetime),':','-'));
% save(filName);
% savefig(filName);
% exportgraphics(gcf,[filName,'.png'],'Resolution',600)
% end



