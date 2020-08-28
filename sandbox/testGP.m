clear
clc
close all

%% generate synthetic flow data
% altitudes
altitudes = 0:10:100;
% altitude length scale
altitudeLenghtScale = 10;
% final time for data generation in minutes
tFinData = 30;
% temporal length scale
temporalLengthScale = 10;
% time step for synthetic data generation
timeStepSynData = 1;
% standard deviation for synthetic data generation
stdDevSynData = 1;
% get the time series object
synFlow = generateSyntheticWindData(altitudes,altitudeLenghtScale,...
    temporalLengthScale,tFinData,stdDevSynData,'timeStep',timeStepSynData);

%% initialize 



%% create gp class object
gp = GP.GaussianProcess(1,'squaredExponential','squaredExponential');
% training points
tpIndices = randperm(numel(X),20);
xTrain = X(tpIndices);
tTrain = T(tpIndices);
yTrain = Y(tpIndices)';

%% optimize hyper-parameters
gp.noiseVariance = 0.001;
XTtrain = [xTrain; tTrain];
iniGuess = [1;1;1;1];
optHyp = gp.findOptSpatioTemporalHyperParams(XTtrain,yTrain,iniGuess);
gp = gp.setOptimumHyperParameters(optHyp);


%% gaussian process regression
xQuery = X(:)';
tQuery = T(:)';
xtQuery = [xQuery; tQuery];

mesh(X,T,Y,'edgecolor','k');
hold on; grid on;
set(gcf,'Position',[-9 133 560 420]);

% number of iterations
noIter = 6;
covMat = gp.makeTotalCovarianceMatrix(XTtrain);
for ii = 1:noIter
    if ii > 1
        delete(scPlot);
    end
    
    % not visited points
    notVisited = setdiff(1:numel(X),tpIndices);
    % visit new point
    newPtIdx = randperm(numel(notVisited),10);
    tpIndices = [tpIndices newPtIdx];
    XTNew     = [X(notVisited(newPtIdx)); T(notVisited(newPtIdx))];
    yNew      = Y(notVisited(newPtIdx))';
    % augment covariance matrix
    covMat = gp.augmentCovarianceMatrix(XTtrain,XTNew,covMat);
    XTtrain = [XTtrain XTNew];
    yTrain  = [yTrain; yNew];
    if ii < 6
        optHyp = gp.findOptSpatioTemporalHyperParams(XTtrain,yTrain,iniGuess);
        gp = gp.setOptimumHyperParameters(optHyp);
    end
    % regression
    [predMean,postVar] = gp.calcPredMeanAndPostVar(covMat,XTtrain,...
        yTrain,xtQuery);
    
    scPlot = scatter3(xQuery,tQuery,predMean,'ro');
    keyboard
    
end





