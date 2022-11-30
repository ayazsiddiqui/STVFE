clear
clc
% close all

cd(fileparts(mfilename('fullpath')));


%% initialize KFGP
rng(1);

% altitudes
altitudes = 0:100:1000;
kfgpTimeStep = 0.5;

% make class object
spKernel = 'squaredExponential';
tKernel  = 'squaredExponential';
kfgp = GP.KalmanFilteredGaussianProcess(spKernel,tKernel,...
    'zeroMean',altitudes,kfgpTimeStep);

kfgp.spatialCovAmp       = 5.1^2;
kfgp.spatialLengthScale  = 220;
kfgp.temporalCovAmp      = 1;
kfgp.temporalLengthScale = 22;
kfgp.noiseVariance       = 0.1;

kfgp.initVals = kfgp.initializeKFGP;
kfgp.spatialCovMat = kfgp.makeSpatialCovarianceMatrix(altitudes);
kfgp.spatialCovMatRoot = kfgp.calcSpatialCovMatRoot;


% guassian process
gp = GP.GaussianProcess(spKernel,tKernel,'zeroMean');

gp.spatialCovAmp       = kfgp.spatialCovAmp;
gp.spatialLengthScale  = kfgp.spatialLengthScale;
gp.temporalCovAmp      = kfgp.temporalCovAmp;
gp.temporalLengthScale = kfgp.temporalLengthScale;
gp.noiseVariance       = kfgp.noiseVariance;

%% generate synthetic flow data
gp2 = GP.GaussianProcess(spKernel,tKernel,'windPowerLaw');
gp2.spatialCovAmp       = kfgp.spatialCovAmp;
gp2.spatialLengthScale  = kfgp.spatialLengthScale;
gp2.temporalCovAmp      = kfgp.temporalCovAmp;
gp2.temporalLengthScale = kfgp.temporalLengthScale;
gp2.noiseVariance       = kfgp.noiseVariance;
% number of altitudes
nAlt = numel(altitudes);
% final time for data generation in minutes
tFinData = 300;
% time step for synthetic data generation
timeStepSynData = 1;
% standard deviation for synthetic data generation
stdDevSynData = 4;
% get the time series object
[synFlow,synAlt] = gp2.generateSyntheticFlowData(altitudes,tFinData,stdDevSynData,...
    'timeStep',timeStepSynData);

%% regression using traditional GP
% algorithm final time
algFinTime = 180;
% sampling time vector
tSamp = 0:kfgp.kfgpTimeStep:algFinTime;
% number of samples
nSamp = numel(tSamp);

% preallocat sampling matrices
xSamp   = NaN(1,nSamp);
flowVal = NaN(1,nSamp);
ySamp   = NaN(nSamp,1);
XTSamp  = NaN(2,nSamp);

nfinAlt = 1*numel(altitudes);
finAlt = linspace(min(altitudes),max(altitudes),nfinAlt);

meanFnVec = gp2.meanFunction(finAlt);

% preallocate matrices for GP
predMeansGP = NaN(nfinAlt,nSamp);
postVarsGP  = NaN(nfinAlt,nSamp);
timeKFGP      = NaN(1,nSamp);
RMSEFitKFGP = NaN(1,nSamp);
stdDevGP    = NaN(nfinAlt,nSamp);
upBoundGP   = NaN(nfinAlt,nSamp);
loBoundGP   = NaN(nfinAlt,nSamp);

% preallocate matrices for KFGP
predMeansKFGP = NaN(nfinAlt,nSamp);
postVarsKFGP  = NaN(nfinAlt,nSamp);
timeGP        = NaN(1,nSamp);
RMSEFitGP     = NaN(1,nSamp);
stdDevKFGP    = NaN(nfinAlt,nSamp);
upBoundKFGP   = NaN(nfinAlt,nSamp);
loBoundKFGP   = NaN(nfinAlt,nSamp);


kfgpToGpFit   = NaN(1,nSamp);
% number of std deviations for bounds calculations
numStdDev = 2;

% make for loop
for ii = 1:nSamp
    % go to xSamp
    if ii == 1
        xSamp(ii) = altitudes(randperm(nAlt,1));
    else
        [~,maxVarIdx] = max(postVarsKFGP(:,ii-1));
        xSamp(ii) = finAlt(maxVarIdx);
    end
    % measure flow at xSamp(ii) at tSamp(ii)
    fData = resample(synFlow,tSamp(ii)*60).Data;
    hData = resample(synAlt,tSamp(ii)*60).Data;
    flowVal(ii) = interp1(hData,fData,xSamp(ii));
    ySamp(ii) =  gp2.meanFunction(xSamp(ii)) - flowVal(ii);
    % augment altitude and height in XTsamp
    XTSamp(:,ii) = [xSamp(ii);tSamp(ii)];
    % recursion
    tic
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
    [muKFGP,sigKFGP] = kfgp.calcPredMeanAndPostVar(finAlt,F_t,sigF_t);  
    % KFGP: store them
    predMeansKFGP(:,ii) = meanFnVec(:) - muKFGP;
    postVarsKFGP(:,ii)  = sigKFGP;
    timeKFGP(ii) = toc;
    % percentage fit to data
    RMSEFitKFGP(ii) = mean(((predMeansKFGP(:,ii) - ...
        interp1(hData,fData,nfinAlt))/mean(meanFnVec)).^2)^0.5; 
    RMSEFitKFGP(ii) = 1 - (norm(predMeansKFGP(:,ii) - fData(:))/max(norm(fData(:)),eps));
    % KFGP: calculate bounds
    stdDevKFGP(:,ii) = postVarsKFGP(:,ii).^0.5;
    % KFGP: upper bounds = mean + x*(standard deviation)
    upBoundKFGP(:,ii) = predMeansKFGP(:,ii) + numStdDev*stdDevKFGP(:,ii);
    % KFGP: lower bounds = mean - x*(standard deviation)
    loBoundKFGP(:,ii) = predMeansKFGP(:,ii) - numStdDev*stdDevKFGP(:,ii);
    tic
    if ii == 1
        % GP: covariance matrix
        covMat = gp.makeTotalCovarianceMatrix(XTSamp(:,ii));        
    else
        % GP: covariance matrix
        covMat = gp.augmentCovarianceMatrix(XTSamp(:,1:ii-1),...
            XTSamp(:,ii),covMat);
    end
    % GP: calculate prediction mean and posterior variance
    [muGP,sigGP] = ...
        gp.calcPredMeanAndPostVar(covMat,XTSamp(:,1:ii),ySamp(1:ii),...
        [finAlt;tSamp(ii)*ones(1,nfinAlt)]);
    % GP: store them
    predMeansGP(:,ii) = meanFnVec(:) - muGP;
    postVarsGP(:,ii)  = sigGP;
    timeGP(ii) = toc;
    % percentage fit to data
    RMSEFitGP(ii) = mean(((predMeansGP(:,ii) - ...
        interp1(hData,fData,nfinAlt))/mean(meanFnVec)).^2)^0.5; 
    RMSEFitGP(ii) = 1 - (norm(predMeansGP(:,ii) - fData(:))/max(norm(fData(:)),eps)); 
    
    kfgpToGpFit(ii) = 1 - (norm(predMeansGP(:,ii) - predMeansKFGP(:,ii))/...
        max(norm(predMeansGP(:,ii)),eps)); 
    % GP: calculate bounds
    stdDevGP(:,ii) = postVarsGP(:,ii).^0.5;
    % GP: upper bounds = mean + x*(standard deviation)
    upBoundGP(:,ii) = predMeansGP(:,ii) + numStdDev*stdDevGP(:,ii);
    % GP: lower bounds = mean - x*(standard deviation)
    loBoundGP(:,ii) = predMeansGP(:,ii) - numStdDev*stdDevGP(:,ii);
    
end

%% plot
cols = [228,26,28
55,126,184
77,175,74
152,78,163
255,127,0
166,86,40
231,41,138]/255;

P.linProp = {...
    cols(1,:),'-';
    cols(2,:),':';
    cols(3,:),'-^';
    cols(4,:),'--';
    cols(5,:),'-o';
    cols(6,:),'-s';
    cols(7,:),'-d'};

tSampPlot = tSamp/60;

% create axis object
nSp = 1;
sbAxis(nSp) = subplot(2,1,nSp);
% set axis properties
hold(sbAxis(nSp),'on');
ylabel(sbAxis(nSp),'Fit [\%]');
fPlots(1) = plot(tSampPlot,100*kfgpToGpFit...
       ,P.linProp{1,2}...
       ,'color','k'...
       ,'MarkerFaceColor',P.linProp{1,1}...
       ,'linewidth',1.2);

nSp = nSp + 1;
sbAxis(nSp) = subplot(2,1,nSp);

fPlots(2) = semilogy(tSampPlot,timeGP...
       ,P.linProp{1,2}...
       ,'color',P.linProp{1,1}...
       ,'MarkerFaceColor',P.linProp{2,1}...
       ,'linewidth',1.2);   
   
hold(sbAxis(nSp),'on');
ylabel(sbAxis(nSp),'Computation time [s]');

fPlots(2) = semilogy(tSampPlot,timeKFGP...
       ,P.linProp{1,2}...
       ,'color',P.linProp{2,1}...
       ,'MarkerFaceColor',P.linProp{4,1}...
       ,'linewidth',1.2);     
ylim([1e-4 1e-1]);

grid(sbAxis(1:end),'on');
set(sbAxis(1:end),'GridLineStyle',':')
xlabel(sbAxis(1:end),'Time [hr]');   
sbAxis(1).YLim = [90 110];

set(sbAxis(1:end),'FontSize',12);
set(gcf,'InnerPosition',1*[-00 -00 560 1.3*420])
sbAxis(1).YLabel.Position(1) = sbAxis(2).YLabel.Position(1);
legend('Batch GP','Kalman filter GP','location','south','orientation','horizontal')

%% export file
saveFile = input('Save file? Options: Enter y or n\n','s');
if strcmpi(saveFile,'y')
filName = strcat('gpAndKfgpComp_',strrep(datestr(datetime),':','-'));
savefig(filName);
exportgraphics(gcf,[filName,'.png'],'Resolution',600)
end


%% convert results to time series and store in strcut
regressionRes(1).predMean  = timeseries(predMeansKFGP,tSamp*60);
regressionRes(1).loBound   = timeseries(loBoundKFGP,tSamp*60);
regressionRes(1).upBound   = timeseries(upBoundKFGP,tSamp*60);
regressionRes(1).dataSamp  = timeseries([xSamp;flowVal],tSamp*60);
regressionRes(1).dataAlts  = timeseries(repmat(finAlt(:),1,nSamp),tSamp*60);
regressionRes(1).legend    = 'KFGP';

regressionRes(2).predMean  = timeseries(predMeansGP,tSamp*60);
regressionRes(2).loBound   = timeseries(loBoundGP,tSamp*60);
regressionRes(2).upBound   = timeseries(upBoundGP,tSamp*60);
regressionRes(2).dataSamp  = timeseries([xSamp;flowVal],tSamp*60);
regressionRes(2).dataAlts  = timeseries(repmat(finAlt(:),1,nSamp),tSamp*60);
regressionRes(2).legend    = 'GP';

%% plot the data
figure
F = animatedPlot(synFlow,synAlt,'plotTimeStep',kfgpTimeStep,...
    'regressionResults',regressionRes,'wait',false,'plotTimeStep',1);

%% video
% % % % video setting
video = VideoWriter('kfgpDemo','Uncompressed AVI');
% % video = VideoWriter('vid_Test1','MPEG-4');
video.FrameRate = 5;
% video.Quality   = 100;
% set(gca,'nextplot','replacechildren');

open(video)
for ii = 1:length(F)
    writeVideo(video, F(ii));
end
close(video)


