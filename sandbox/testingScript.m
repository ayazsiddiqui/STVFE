clear
clc

%% make instance of the gpkf
FastGPKF = GPKF(1,"exponential","upperConfidenceBound");

% set hyperparameters
FastGPKF.p_temporalLenghtScale = 10;    % minutes
FastGPKF.p_spatialCovarianceAmp = 1;    % ul
FastGPKF.p_noiseVariance = 0.01;        %
FastGPKF.p_spatialLenghtScale = 200;    % meters

% set gpkf parameters
FastGPKF.p_gpfkTimeStep = 0.5;
FastGPKF.p_xMeasure = 100:100:500;

% set acquisition function parameters
FastGPKF.p_exploitationConstant = 0;
FastGPKF.p_explorationConstant = 2;

% set initial values
FastGPKF = FastGPKF.m_CalcSpatialCovMat;
FastGPKF = FastGPKF.m_CalcSpatialCovRoot;
FastGPKF = FastGPKF.gpfkInitialize;

%% test methods
% dummy values
x1 = [100 200;10 20];
y1 = [5;10];
x2 = 2.*x1;
y2 = 2.*y1;

% test m_CalcSpatialCovariance
test_1 = FastGPKF.m_CalcSpatialCovariance(x1(1,1),x2(1,1));

% test m_buildSpatialCovMat
test_2 = FastGPKF.m_buildSpatialCovMat(x1(1,:));

% test m_temporalCovariance
test_3 = FastGPKF.m_temporalCovariance(x1(2,1),x2(2,1));

% test m_calcTotCovariance
test_4 = FastGPKF.m_calcTotCovariance(x1(:,1),x2(:,1));

% test m_buildCovMatAndMeanVec
[test_5.covMat,test_5.meanVec] = FastGPKF.m_buildCovMatAndMeanVec(x1);

% test m_calcMarginalLikelihood
test_6 = FastGPKF.m_calcMarginalLikelihood(x1,y1);

% test m_gpkfInitialize
test_7 = FastGPKF.m_gpkfInitialize();

% test m_gpkfKalmanEstimation
sk_k = FastGPKF.p_gpfkInitVals.s0;
ck_k = FastGPKF.p_gpfkInitVals.sig0Mat;
Mk = FastGPKF.p_xMeasure(2);

[F_t,sigF_t,skp1_kp1,ckp1_kp1,Ik] = FastGPKF.m_gpkfKalmanEstimation...
    (sk_k,ck_k,Mk,y1(1));

% test m_gpkfRegression
xPredict = linspace(FastGPKF.p_xMeasure(1),FastGPKF.p_xMeasure(end),...
    2*numel(FastGPKF.p_xMeasure)-1);
[predMean,postVar] = FastGPKF.m_gpkfRegression(xPredict,F_t,sigF_t);

% test m_predictionGPKF
Mk2 = FastGPKF.p_xMeasure(randperm(numel(FastGPKF.p_xMeasure),3));
op = FastGPKF.m_predictionGPKF(sk_k,ck_k,Mk,y1(1),Mk2,3);

% test m_gpkfMPC_bruteForce
uAllowable = -100:100:100;
testBruteForce = FastGPKF.m_gpkfMPC_bruteForce(sk_k,ck_k,Mk,y1(1),...
    uAllowable,4);

% test m_gpkfMPC_fmincon
testFmincon = FastGPKF.m_gpkfMPC_fmincon(sk_k,ck_k,Mk,y1(1),...
    uAllowable,4,5);




