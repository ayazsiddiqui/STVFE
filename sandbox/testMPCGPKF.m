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
timeStep = 0.05*10;
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
% % % prediction horizon
predHorz = 4;
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

% initCons = GpkfInitialize(xDomain,...
%     optHyperParams(end-1),timeStep,'approximationOrder',Nn);

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
predMean = NaN(size(xPredict,2),predHorz,noIter);
postVar =  NaN(size(xPredict,2),predHorz,noIter);
stdDev =  NaN(size(xPredict,2),predHorz,noIter);
upperBound = NaN(size(xPredict,2),predHorz,noIter);
lowerBound = NaN(size(xPredict,2),predHorz,noIter);
pointsVisited = NaN(nVisit,noIter);
fValAtPt = NaN(noIter,nVisit);
GPKFcompTime = NaN(1,noIter);
GPKFfit = NaN(1,noIter);

for ii = 1:noIter
    % % % visit said points
    if ii == 1
        visitIdx = sort(randperm(size(xMeasure,2),nVisit));
    else
        [~,visitIdx] = max(sum(postVar(:,:,ii-1),2));
    end
    % % % extract visited values from xMeasure
    Mk = xMeasure(:,visitIdx);
    % % % extract wind speed at visited values
    yk = windSpeedOut(visitIdx,ii);
    % % % stepwise update of kalman state estimate and covariance
    tic
    [F_t,sigF_t,skp1_kp1,ckp1_kp1] = ...
        gpkf.gpkfKalmanEstimation(xMeasure,sk_k,ck_k,Mk,yk,...
        Ks_12,initCons.Amat,initCons.Qmat,initCons.Hmat,optHyperParams(end),...
        predHorz);
    % % % regression over a finer domain
    for jj = 1:predHorz
        [predMean(:,jj,ii),postVar(:,jj,ii)] = gpkf.gpkfRegression(xDomain,xPredict,...
            F_t(:,jj),sigF_t(:,:,jj),Ks,optHyperParams);
        % % % remove real or imaginary parts lower than eps
        stdDev(:,jj,ii) = sqrt(gpkf.removeEPS(postVar(:,jj,ii),5));
        % % % upper bounds = mean + x*(standard deviation)
        upperBound(:,jj,ii) = predMean(:,jj,ii) + 1*stdDev(:,jj,ii);
        % % %lower bounds = mean + x*(standard deviation)
        lowerBound(:,jj,ii) = predMean(:,jj,ii) - 1*stdDev(:,jj,ii);
    end
    % % % store points visited at the respective function value
    pointsVisited(:,ii) = Mk(:);
    fValAtPt(ii,:) = yk(:);
    % % % update previous step information
    sk_k = skp1_kp1(:,1);
    ck_k = ckp1_kp1(:,:,1);
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

for ii = 1:noTimeSteps-predHorz+1
    
    for jj = 1:predHorz
        subplot(2,2,jj)
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
    end
    
    % % plot true wind
    for jj = 1:predHorz
        subplot(2,2,jj)
        plTrueWind = plot(windSpeedOut(:,ii+jj-1),heights,'k','linewidth',lwd);
        % % plot measured wind value
        plfVals = plot(fValAtPt(ii,:),pointsVisited(:,ii),'mo',...
            'markerfacecolor','m','linewidth',lwd);
        % % plot GPKF mean and bounds
        plPredMean = plot(predMean(:,jj,ii),xPredict,'-x','linewidth',lwd,...
            'color',1/255*[228,26,28]);
        plLowerBds = plot(lowerBound(:,jj,ii),xPredict,'--','linewidth',lwd,...
            'color',1/255*[254,178,76]);
        plUpperBds = plot(upperBound(:,jj,ii),xPredict,'--','linewidth',lwd,...
            'color',1/255*[254,178,76]);
        % % legend
        legend([plTrueWind,plPredMean,plLowerBds],...
            'True func','GPKF $\mu$','GPKF bounds');
        % % title
        txt1 = sprintf('$l_{t} / \\tau$ = %0.2f,',timeScale/timeStep);
        txt = sprintf(' Time = %0.2f min',tVec(ii+jj-1));
        txt = strcat(txt1,txt);
        title(txt);
    end
    
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
% % % % video setting
% video = VideoWriter(fName,'Motion JPEG AVI');
% % % video = VideoWriter('vid_Test1','MPEG-4');
% video.FrameRate = 10;
% set(gca,'nextplot','replacechildren');
%
% open(video)
% for ii = 1:length(F)
%     writeVideo(video, F(ii));
% end
% close(video)


