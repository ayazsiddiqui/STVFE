clear
clc
format compact
close all

%% change working directory to current folder
cd(fileparts(mfilename('fullpath')));

%% generate wind using colored noise
rng(56);
% environment
hMax = 1500;
hMin = 100;
heights = hMin:100:hMax;
heights = heights(:);
meanFlow = 10;
noTP = numel(heights);
% time in minutes
timeStep = 0.05*2;
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
windSpeedOut = 1*(0+genWindv2(heights,heightScale,tVec,timeScale,stdDev));
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
gpkf = GPKF(1,'exponential','upperConfidenceBound');
% set values of hyper parameters
noiseVar = 0.01;
hyperParams = [1 heightScale timeScale noiseVar]';
% introduce some error in the hyper parameters
misCalc = (0.5-0).*rand(4,1);
optHyperParams = 1.*hyperParams;
% quantify the error
hyperError = optHyperParams./hyperParams;

%% set up MPC
% % % prediction horizon
predHorz = 1;
% % % time step for mpc trigger
mpcInterval = timeScale/5;
% % % allowable control inputs
uAllowable = -1000:100:1000;
% % % exploration constant: used as (2^(c))*(posterior variance)
exploreConstant = 0;
% % % exploitation constant: used as c*(prediction mean)
exploitConstant = 0;

%% do the actual gaussian process kalman filtering
% % % make domain vector
xDomain = heights(:)';
% % % set the measurable domain equal to the entire domain
xMeasure = xDomain;
% % % make a finer domain over which predictions are made
xPredict = linspace(heights(1),heights(end),1*numel(heights));
% % % order of approixation for SE kernel
Nn = 2;
% form the initialization matrices
initCons = gpkf.GpkfInitialize(xDomain,...
    optHyperParams(end-1),timeStep,'approximationOrder',Nn);
% form initial matrices for MPC
initConsMPC = gpkf.GpkfInitialize(xDomain,...
    optHyperParams(end-1),mpcInterval,'approximationOrder',Nn);

% % % initialize parameters
Ks = gpkf.buildSpatialCovMat(xMeasure,optHyperParams(1),optHyperParams(2));

% Ks_12 = chol(Ks,'upper');
% Ks_12 = Ks_12 + triu(Ks_12,1)';
Ks_12 = sqrtm(Ks);

ck_k = initCons.sig0Mat;
sk_k = initCons.s0;

% % % number of iterations
noIter = noTimeSteps;
% % % preallocate matrices
predMean = NaN(size(xPredict,2),noIter);
postVar =  NaN(size(xPredict,2),noIter);
acquiFun = NaN(size(xPredict,2),noIter);
stdDev =  NaN(size(xPredict,2),noIter);
upperBound = NaN(size(xPredict,2),noIter);
lowerBound = NaN(size(xPredict,2),noIter);
pointsVisited = NaN(1,noIter);
fValAtPt = NaN(noIter,1);
GPKFcompTime = NaN(1,noIter);
GPKFfit = NaN(1,noIter);
optCtrlSequence = NaN(ceil(tVec(end)/mpcInterval),predHorz);
% start off mpc counter
mpcCount = 1;

for ii = 1:noIter
    % % % visit said points
    if ii == 1 || ~exist('mpcRes','var')
        visitIdx = randperm(size(xMeasure,2),1);
        Mk = xMeasure(:,visitIdx);
    else
        [~,maxIdx] =  max(postVar(:,ii-1));
        Mk2 = xMeasure(:,maxIdx);
        Mk = mpcRes.optStateTrajectory(2);
    end
    % % % extract wind speed at visited values
    yk = windSpeedOut((Mk == xMeasure)',ii);
    % % % gpkf MPC
    if mod(tVec(ii),mpcInterval) == 0
        mpcRes = gpkf.gpkfMPC(xMeasure,sk_k,ck_k,Mk,yk,Ks_12,...
            initConsMPC.Amat,initConsMPC.Qmat,initConsMPC.Hmat,xPredict,...
            Ks,optHyperParams,uAllowable,predHorz,...
            'explorationConstant',exploreConstant,...
            'exploitationConstant',exploitConstant);
        % % % store optimal control sequence)
        optCtrlSequence(mpcCount,:) = mpcRes.optCtrlSeq;
        mpcCount = mpcCount + 1;
    end
    % % % kalman state estimation
    [F_t,sigF_t,skp1_kp1,ckp1_kp1] = ...
        gpkf.gpkfKalmanEstimation(xMeasure,sk_k,ck_k,Mk,yk,...
        Ks_12,initCons.Amat,initCons.Qmat,initCons.Hmat,...
        optHyperParams(end));
    % % % regression over a finer domain
    [predMean(:,ii),postVar(:,ii)] = gpkf.gpkfRegression(xDomain,xPredict,...
        F_t,sigF_t,Ks,optHyperParams);
    % % % remove real or imaginary parts lower than eps
    stdDev(:,ii) = postVar(:,ii).^0.5;
    % % % upper bounds = mean + x*(standard deviation)
    upperBound(:,ii) = predMean(:,ii) + 1*stdDev(:,ii);
    % % % lower bounds = mean + x*(standard deviation)
    lowerBound(:,ii) = predMean(:,ii) - 1*stdDev(:,ii);
    % % % calculate acquisition function value
    acquiFun(:,ii) = gpkf.calcAcquisitionFun(predMean(:,ii),postVar(:,ii),...
        max(predMean(:,ii)),...
        'explorationConstant',exploreConstant,...
        'exploitationConstant',exploitConstant);
    
    
    
    % % % update previous step information
    sk_k = skp1_kp1;
    ck_k = ckp1_kp1;
    
    %     op = gpkf.predictionGPKF(xMeasure,sk_k,ck_k,Mk,yk,...
    %                 Ks_12,initCons.Amat,initCons.Qmat,initCons.Hmat,...
    %                 xPredict,Ks,optHyperParams,predHorz);
    
    % % % store points visited at the respective function value
    pointsVisited(:,ii) = Mk(:);
    fValAtPt(ii,:) = yk(:);
    
end

%% plot data
% keyboard
% % % linewidths
lwd = 1;
% % % set find plot limits
lB = min([lowerBound(:);windSpeedOut(:)]);
uB = max([upperBound(:);windSpeedOut(:)]);
plotRes = 1;

figure(1)
set(gcf,'Units','normalized','position',[1 0.0889 0.6667 0.8398]);
F = struct('cdata',uint8(zeros(840,1680,3)),'colormap',[]);

for ii = 1:noTimeSteps
    
    %     for jj = 1:predHorz
    %         subplot(2,2,jj)
    if ii == 1
        hold on
        grid on
        xlabel('Wind speed (m/s)');
        ylabel('Altitude (m)');
        %         xlim([lB-mod(lB,plotRes),uB-mod(uB,plotRes)+plotRes])
        xlim(0+[-4 4])
        ylim([hMin hMax]);
    else
        delete(findall(gcf,'type','annotation'));
        h = findall(gca,'type','line','linestyle','-','-or',...
            'linestyle','--','-or','linestyle','-x','-or',...
            'linestyle','-.', '-or','color','m');
        delete(h);
        
    end
    %     end
    
    % % plot true wind
    %     for jj = 1:predHorz
    %         subplot(2,2,jj)
    plTrueWind = plot(windSpeedOut(:,ii),heights,'k','linewidth',lwd);
    % % plot measured wind value
    plfVals = plot(fValAtPt(ii,:),pointsVisited(:,ii),'mo',...
        'markerfacecolor','m','linewidth',lwd);
    % % plot GPKF mean and bounds
    plPredMean = plot(predMean(:,ii),xPredict,'-x','linewidth',lwd,...
        'color',1/255*[228,26,28]);
    plLowerBds = plot(lowerBound(:,ii),xPredict,'--','linewidth',lwd,...
        'color',1/255*[254,178,76]);
    plUpperBds = plot(upperBound(:,ii),xPredict,'--','linewidth',lwd,...
        'color',1/255*[254,178,76]);
    % % legend
    legend([plTrueWind,plPredMean,plLowerBds],...
        'True func','GPKF $\mu$','GPKF bounds');
    % % title
    txt1 = sprintf('$l_{t} / \\tau$ = %0.2f,',timeScale/timeStep);
    txt = sprintf(' Time = %0.2f min',tVec(ii));
    txt = strcat(txt1,txt);
    title(txt);
    %     end
    ff = getframe(gcf);
    F(ii).cdata = ff.cdata;
    
end

%% other plots
% % % % computational time plot
figure(2)
x = gcf;
set(gcf,'position',x.Position.*[1 0 1 1])
hold on
grid on
pGPKFC = plot(1:noIter,GPKFcompTime,'linewidth',lwd,...
    'color',1/255*[228,26,28]);
ylim([-0.005 Inf])
xlabel('Time step number')
ylabel('Computational time (sec)')
legend(pGPKFC,{'GPKF'},'location','best')

figure(3)
x = gcf;
set(gcf,'position',x.Position.*[1 0 1 1])
hold on
grid on
pGPKFFit = plot(tVec,GPKFfit,'linewidth',lwd,...
    'color',1/255*[228,26,28]);
xlabel('Time (min)')
ylabel('Fit (\%)')
ylim([0 100])
legend(pGPKFFit,{'GPKF'},'location','best')

set(findobj('-property','FontSize'),'FontSize',12)

%% save data to output folder
[status, msg, msgID] = mkdir(pwd,'outputs');
fName = [pwd,'\outputs\',strrep(datestr(datetime),':','_')];

% delete([pwd,'\outputs\*.mat'])
% delete([pwd,'\outputs\*.avi'])

% save(fName)

%% video
% % % video setting
video = VideoWriter(fName,'Motion JPEG AVI');
% % video = VideoWriter('vid_Test1','MPEG-4');
video.FrameRate = 10;
set(gca,'nextplot','replacechildren');

open(video)
for ii = 1:length(F)
    writeVideo(video, F(ii));
end
close(video)


