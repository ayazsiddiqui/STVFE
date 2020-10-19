clear all 
close all
clc

load('SimOutRGP.mat')

basisDesignSpace = linspace(0,20,size_grid)*(pi/180);

thetaLib = [2 4 3 1 0 15 5 7 6 8 12 8 8 8]'*(pi/180);
PerfmIndxLib = [0.3 0.22 0.35 0.18 0.4 0.33 0.42 0.24 0.26 0.19 0.25 0.25 0.22 0.22]';
% sigmaSq = 0.05;        % Approximate variance level due to noise
sigmaSq = 1;        % Approximate variance level due to noise
hypPara = 0.01;
predVarOfflineSimVec = ones(size_grid,2);

for ii = 1:length(basisDesignSpace)
    for jj = 1:length(basisDesignSpace)
        basisCov(ii,jj) = covFuncEval(basisDesignSpace(ii),basisDesignSpace(jj),hypPara,sigmaSq); 
    end 
end 

% Initialize covaraince and mean function
invbasisCov = ones(size_grid,size_grid)/(basisCov);
covUpdateF_prev = basisCov;
meanFuncF_prev = 0.0*ones(length(basisDesignSpace),1);
meanFuncBasis = 0.0*ones(length(basisDesignSpace),1);

%-------------------------------------------------------------------------
offlineSimLength = 10;      % Length time considered for this offline simulation
offlineSimEndSample = 10/sample_time + 1;

% Pick off a sequence of points from the simulation data
sampleAdjustment = 60/sample_time;            % This adjustment is due to time during the simulation where I allow the model to settle
thetaLibRad = thetaMeas(sampleAdjustment:sampleAdjustment+offlineSimEndSample,1)*(pi/180);
PerfmIndxLib = instPerfmIndex(sampleAdjustment:sampleAdjustment+offlineSimEndSample,1);

filtResponseValue0 = instPerfmIndex(sampleAdjustment-1,1);

for jj = 1:offlineSimEndSample
    
    desParaCur = thetaLibRad(jj,1);
    responseCurr = PerfmIndxLib(jj,1);
    
    for ii = 1:size_grid
        covCurPara(1,ii) = covFuncEval(desParaCur,basisDesignSpace(ii),hypPara,sigmaSq);
    end 
    covCurrPoint = covFuncEval(desParaCur,desParaCur,hypPara,sigmaSq);
    
    % Distance calculation to linearly interpolate between nearest 2 points
    % to the current point
    distCur2Grid = sort(abs(basisDesignSpace-desParaCur));   
    firstDistIndex = find(abs(basisDesignSpace-desParaCur) == distCur2Grid(1));
    secDistIndex = find(abs(basisDesignSpace-desParaCur) == distCur2Grid(2));
    if length(secDistIndex) == 2
        % If the current point corresponds to a grid point, use the mean at
        % that point, not an interpolation between 2 points
        trueIndex = round(mean(secDistIndex));
        meanCurPara = meanFuncF_prev(trueIndex);
    else 
        CurParaInterp = [basisDesignSpace(firstDistIndex) basisDesignSpace(secDistIndex)];
        meanCurParaInterp = [meanFuncF_prev(firstDistIndex) meanFuncF_prev(secDistIndex)];
        % Acutal function call to linearly interpolate between closest points
        meanCurPara = interp1(CurParaInterp,meanCurParaInterp,desParaCur);
    end 
        
    % Inference part of the algorithm
    covUpdateVec = covCurPara/(basisCov + 1E-6*eye(size_grid,size_grid));   
    covUpdateVecTest = covCurPara*invbasisCov; 
       
    meanFuncPupdate = 0.0 + covUpdateVec*(meanFuncF_prev - meanFuncBasis);
%     meanFuncPupdate = covUpdateVec*meanFuncF_prev;
%     meanFuncPupdate = meanFuncBasis(firstDistIndex) + covUpdateVec*(meanFuncF_prev - meanFuncBasis);
    covUpdatedP = covCurrPoint + covUpdateVec*(covUpdateF_prev - basisCov)*covUpdateVec';
    
    meanFuncPupdateCheck = covUpdateVec*meanFuncF_prev;
    meanFuncPupdateCheckError(jj,1) = meanFuncPupdateCheck - meanFuncPupdate;
    
    covPredictionCheck(jj,1) = covUpdateVec*covUpdateF_prev*covUpdateVec';
    covPredictionCheckError(jj,1) = covPredictionCheck(jj,1) - covUpdatedP;
  
    
    keepCovUpdatedP(jj,1) = covUpdatedP;
    keepCovUpdatedPterm(jj,1) = covUpdateVec*(covUpdateF_prev - basisCov)*covUpdateVec';
    
    tauFiltResponse = 1;
    if jj == 1
        filtResponseValue(jj,1) = (sample_time*responseCurr + tauFiltResponse*filtResponseValue0)/(sample_time + tauFiltResponse);
    else 
        filtResponseValue(jj,1) = (sample_time*responseCurr + tauFiltResponse*filtResponseValue(jj-1))/(sample_time + tauFiltResponse);
    end 
    
    errorFilt(jj,1) = filtResponseValue(jj,1) - responseCurr;
    
    % Update portion of the algorithm
%         kalGainMatrix = covUpdateF_prev*covUpdateVec'*(covUpdatedP + sigmaSq)^-1;
        kalGainMatrix = covUpdateF_prev*covUpdateVec'*(covUpdatedP + sigmaSq)^-1;
        meanFuncFupdate = meanFuncF_prev + kalGainMatrix*(responseCurr - meanFuncPupdate);
        covUpdatedF = covUpdateF_prev - kalGainMatrix*covUpdateVec*covUpdateF_prev;   
        
%         covUpdatedFcheck = covUpdateF_prev - ... 
%             covUpdateVec*covUpdateF_prev*(covUpdatedP+sigmaSq)^(-1)*covUpdateF_prev*covUpdateVec';
        
        covUpdatedFcheck = covUpdateF_prev - ... 
            covUpdateF_prev*covUpdateVec'*(covUpdatedP+sigmaSq)^(-1)*covUpdateVec*covUpdateF_prev;

        maxdiff(jj,1) = max(max(abs(covUpdatedF-covUpdatedFcheck)));
        
%         plot(basisDesignSpace,meanFuncFupdate)
%         grid on
%         hold on
                    
    % Calculation of prediction variance based on updated covariance matrix
    predVarOfflineSim(:,jj) = diag(covUpdatedF);  
    
    % Keep current vector/matrix values for future iterations
    meanFuncF_prev = meanFuncFupdate;
    covUpdateF_prev = covUpdatedF;
    
end 

%--------------------------------------------------------------------------
% Calculation of true predictive mean and variance from traditional GP
% modelling

for ii = 1:offlineSimEndSample
    for jj = 1:offlineSimEndSample
        trueCovMatrix(ii,jj) = covFuncEval(thetaLibRad(ii),thetaLibRad(jj),hypPara,sigmaSq); 
    end 
end 

for ii = 1:size_grid
    for jj = 1:offlineSimEndSample
        covRowVectorTrue(jj,ii) = covFuncEval(basisDesignSpace(ii),thetaLibRad(jj),hypPara,sigmaSq);
    end 
end 


noiseMat = sigmaSq*eye(offlineSimEndSample);
for kk = 1:size_grid
%     trueMeanFunc(kk,1) = covRowVectorTrue(:,kk)'*inv(trueCovMatrix)*PerfmIndxLib(1:offlineSimEndSample);
%     trueCovFunc(kk,1) = 1 - covRowVectorTrue(:,kk)'*inv(trueCovMatrix)*covRowVectorTrue(:,kk);
    trueMeanFunc(kk,1) = (covRowVectorTrue(:,kk)'/(trueCovMatrix + noiseMat))*PerfmIndxLib(1:offlineSimEndSample);
    trueCovFunc(kk,1) = sigmaSq - (covRowVectorTrue(:,kk)'/(trueCovMatrix + noiseMat))*covRowVectorTrue(:,kk);
end 
%--------------------------------------------------------------------------
% Generation of variance window for GP 
meanFuncFUpperLimit = meanFuncFupdate + predVarOfflineSim(:,end);
meanFuncFLowerLimit = meanFuncFupdate - predVarOfflineSim(:,end);
%--------------------------------------------------------------------------
% Debug plotting of mean function
figure
hold on
plot(basisDesignSpace*(180/pi),meanFuncFupdate(:,end))
scatter(thetaLibRad*(180/pi),PerfmIndxLib,'r','filled')
plot(basisDesignSpace*(180/pi),trueMeanFunc,'--k','LineWidth',2)
plot(basisDesignSpace*(180/pi),meanFuncFUpperLimit,'-g','LineWidth',2)
plot(basisDesignSpace*(180/pi),meanFuncFLowerLimit,'-g','LineWidth',2)
xlabel('$\theta$')
ylabel('$J_{inst}$')
grid on
legend('RGP','data','GP')
% legend('RGP','data','GP','$RGP_{UL}$','$RGP_{LL}$')
hold off

% Debug plotting of covariance function
figure
hold on
plot(basisDesignSpace*(180/pi),trueCovFunc,'-k','LineWidth',2)
plot(basisDesignSpace*(180/pi),predVarOfflineSim(:,end),'--r','LineWidth',2)
xlabel('$\theta$')
ylabel('$\sigma ^2 (\theta)$')
legend('GP','RGP')
grid on
hold off

% % Debug plotting of covariance function
% figure
% hold on
% % plot(basisDesignSpace*(180/pi),predVarOfflineSim(:,end))
% plot(basisDesignSpace*(180/pi),trueCovFunc,'-k','LineWidth',2)
% plot(basisDesignSpace*(180/pi),predVarOfflineSimDebug(:,end),'-r','LineWidth',2)
% xlabel('$\theta$')
% ylabel('$\sigma ^2 (\theta)$')
% legend('RGP', 'GP','$RGP"$')
% grid on
% hold off

% figure
% hold on
% plot(basisDesignSpace*(180/pi),predVarOfflineSim(:,end)*sigmaSq + 1)
% plot(basisDesignSpace*(180/pi),trueCovFunc,'-k','LineWidth',2)
% xlabel('$\theta$')
% ylabel('$\sigma ^2 (\theta)$')
% legend('RGP', 'GP')
% grid on
% hold off   
%--------------------------------------------------------------------------    
function [covFuncValue] = covFuncEval(designPt1,designPt2,theta,noiseVar)
    covFuncValue = noiseVar*exp(-1/(2*theta^2)*(abs(designPt1-designPt2))^2);
end 
    