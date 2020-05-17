clear
clc
format compact
close all

%% change working directory to current folder
cd(fileparts(mfilename('fullpath')));

%% generate wind using colored noise
rng(2);

% environment
hMax = 1500;
hMin = 100;
heights = hMin:100:hMax;
heights = heights(:);
meanFlow = 10;
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
windSpeedOut = genWindv2(heights,heightScale,tVec,timeScale,stdDev);
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
gpkf = GPKF(1,'exponential');
% set values of hyper parameters
noiseVar = 0.01;
hyperParams = [1 heightScale timeScale noiseVar]';
% introduce some error in the hyper parameters
misCalc = (0.5-0).*rand(4,1);
optHyperParams = 1.*hyperParams;
% quantify the error
hyperError = optHyperParams./hyperParams;


%% do the actual gaussian process kalman filtering
% % % make domain vector
xDomain = heights(:)';
% % % set the measurable domain equal to the entire domain
xMeasure = xDomain;
% % % make a finer domain over which predictions are made
xPredict = linspace(heights(1),heights(end),1*numel(heights));
% % % order of approixation for SE kernel
Nn = 6;
% form the initialization matrices
initCons = gpkf.GpkfInitialize(xDomain,...
    optHyperParams(end-1),timeStep,'approximationOrder',Nn);
% % % set number of points visited per step
nVisit = 1;
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
stdDev =  NaN(size(xPredict,2),noIter);
upperBound = NaN(size(xPredict,2),noIter);
lowerBound = NaN(size(xPredict,2),noIter);
pointsVisited = NaN(nVisit,noIter);
fValAtPt = NaN(noIter,nVisit);
GPKFcompTime = NaN(1,noIter);
GPKFfit = NaN(1,noIter);

for ii = 1:noIter
    % % % visit said points
    if ii == 1
        visitIdx = sort(randperm(size(xMeasure,2),nVisit));
    else
        [~,sortIdx] =  max(postVar(:,ii-1));
        [~,visitIdx] = min((xMeasure - xPredict(sortIdx)).^2);
    end
    % % % extract visited values from xMeasure
    Mk = xMeasure(:,visitIdx);
    % % % extract wind speed at visited values
    yk = windSpeedOut(visitIdx,ii);
    % % % stepwise update of kalman state estimate and covariance
    tic
    [F_t,sigF_t,skp1_kp1,ckp1_kp1] = ...
        gpkf.gpkfKalmanEstimation(xMeasure,sk_k,ck_k,Mk,yk,...
        Ks_12,initCons.Amat,initCons.Qmat,initCons.Hmat,optHyperParams(end));
    % % % regression over a finer domain
    [predMean(:,ii),postVar(:,ii)] = gpkf.gpkfRegression(xDomain,xPredict,...
        F_t,sigF_t,Ks,optHyperParams);
    GPKFcompTime(ii) = toc;
    % % % percentage fit
    GPKFfit(ii) = 100*(1 - (norm(predMean(:,ii)-windSpeedOut(:,ii))/...
        max([norm(windSpeedOut(:,ii)),eps])));
    % % % remove real or imaginary parts lower than eps
    stdDev(:,ii) = sqrt(gpkf.removeEPS(postVar(:,ii),5));
    % % % upper bounds = mean + x*(standard deviation)
    upperBound(:,ii) = predMean(:,ii) + 1*stdDev(:,ii);
    % % %lower bounds = mean + x*(standard deviation)
    lowerBound(:,ii) = predMean(:,ii) - 1*stdDev(:,ii);
    % % % store points visited at the respective function value
    pointsVisited(:,ii) = Mk(:);
    fValAtPt(ii,:) = yk(:);
    % % % update previous step information
    sk_k = skp1_kp1;
    ck_k = ckp1_kp1;
end

%% tradtional GP estimation
% % % preallocate matrices
GPpredMean = NaN(size(xPredict,2),noIter);
GPpostVar =  NaN(size(xPredict,2),noIter);
GPstdDev =  NaN(size(xPredict,2),noIter);
GPupperBound = NaN(size(xPredict,2),noIter);
GPlowerBound = NaN(size(xPredict,2),noIter);
GPcompTime = NaN(1,noIter);
GPfit = NaN(1,noIter);

for ii = 1:noIter
    % % % all points visited and y at values upto step ii
    xVisited = [pointsVisited(:,1:ii);tVec(1,1:ii)];
    yVisited = fValAtPt(1:ii,:);
    % % % time at which prediction is desired
    tPredict = tVec(ii);
    % % % traditional GP regression
    tic
    [GPpredMean(:,ii),GPpostVar(:,ii)] = gpkf.traditionalGpRegression(xVisited,...
        yVisited,xPredict,tPredict,optHyperParams);
    GPcompTime(ii) = toc;
    % % % percentage fit
    GPfit(ii) = 100*(1 - (norm(GPpredMean(:,ii)-windSpeedOut(:,ii))/...
        max([norm(windSpeedOut(:,ii)),eps])));
    % % % remove real or imaginary parts lower than eps
    GPstdDev(:,ii) = sqrt(gpkf.removeEPS(GPpostVar(:,ii),5));
    % % % upper bounds = mean + x*(standard deviation)
    GPupperBound(:,ii) = GPpredMean(:,ii) + 1*GPstdDev(:,ii);
    % % % lower bounds = mean + x*(standard deviation)
    GPlowerBound(:,ii) = GPpredMean(:,ii) - 1*GPstdDev(:,ii);
    % % % percentage completion
    GPtxt = fprintf('Traditional GP Percentage completion = %0.2f%% \n',...
        100*ii/noIter);
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
set(gcf,'position',[1268 466 1.15*[560 420]]);
F = struct('cdata',uint8(zeros(840,1680,3)),'colormap',[]);

for ii = 1:noTimeSteps
    
    if ii == 1
        hold on
        grid on
        xlabel('Wind speed (m/s)');
        ylabel('Altitude (m)');
%         xlim([lB-mod(lB,plotRes),uB-mod(uB,plotRes)+plotRes])
        xlim([-4 4])
        ylim([hMin hMax]);
    else
        delete(findall(gcf,'type','annotation'));
        h = findall(gca,'type','line','linestyle','-','-or',...
            'linestyle','--','-or','linestyle','-x','-or',...
           'linestyle','-.', '-or','color','m');
        delete(h);
        
    end
    
    % % plot true wind
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
    % % plot GP mean and bounds
    plGPPredMean = plot(GPpredMean(:,ii),xPredict,'-x','linewidth',lwd,...
        'color',1/255*[55,126,184]);
    plGPLowerBds = plot(GPlowerBound(:,ii),xPredict,'-.','linewidth',lwd,...
        'color',1/255*[158,188,218]);
    plGPUpperBds = plot(GPupperBound(:,ii),xPredict,'-.','linewidth',lwd,...
        'color',1/255*[158,188,218]);
    
    legend([plTrueWind,plPredMean,plLowerBds,plGPPredMean,plGPLowerBds],...
        'True func','GPKF $\mu$','GPKF bounds','GP $\mu$','GP bounds')
    txt = sprintf(['$\\frac{Time ~scale}{Time ~step}$ = %0.2f,'...
        'Time = %0.2f min'],...
        timeScale/timeStep,tVec(ii));
    title(txt);
    
    ff = getframe(gcf);
    F(ii).cdata = ff.cdata;
    
end

% % % % computational time plot
figure(2)
x = gcf;
set(gcf,'position',x.Position.*[1 0 1 2])
subplot(2,1,1)
pGPC = plot(1:noIter,GPcompTime,'linewidth',lwd,'color',...
    1/255*[55,126,184]);
hold on
grid on
pGPKFC = plot(1:noIter,GPKFcompTime,'linewidth',lwd,...
    'color',1/255*[228,26,28]);
ylim([-0.01 Inf])
xlabel('Time step number')
ylabel('Computational time (sec)')
legend([pGPC,pGPKFC],{'GP','GPKF'},'location','best')
title(sprintf('Time step = %0.3f sec',timeStep))

subplot(2,1,2)
pGPFit = plot(1:noIter,GPfit,'linewidth',lwd,'color',...
    1/255*[55,126,184]);
hold on
grid on
pGPKFFit = plot(1:noIter,GPKFfit,'linewidth',lwd,...
    'color',1/255*[228,26,28]);
xlabel('Time step number')
ylabel('Fit (\%)')
legend([pGPFit,pGPKFFit],{'GP','GPKF'},'location','best')


%% save data to output folder
% [status, msg, msgID] = mkdir(pwd,'outputs');
% fName = [pwd,'\outputs\',strrep(datestr(datetime),':','_')];
% 
% % delete([pwd,'\outputs\*.mat'])
% % delete([pwd,'\outputs\*.avi'])
% 
% save(fName)
% 
% %% video
% % % % % video setting
% video = VideoWriter(fName,'Motion JPEG AVI');
% % % video = VideoWriter('vid_Test1','MPEG-4');
% video.FrameRate = 3;
% set(gca,'nextplot','replacechildren');
% 
% open(video)
% for ii = 1:length(F)
%     writeVideo(video, F(ii));
% end
% close(video)


