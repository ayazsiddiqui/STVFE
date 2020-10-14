function [predVar,optBasisVec,J_hatStar,covUpdatedF,... 
    meanFuncFupdate]  = RecursiveGPcoloredsnippet_Ayaz(betaBarPrev,AzLimPrev,perfmIndexVal,basisParaGrid,size_grid,... 
    tauFiltColor,sample_time,sigmaSq,hypParaMat,sigmaSqKern,basisCov,invBasisCov,meanFuncIC,...
    covUpdateF_prev,meanFuncF_prev)

predVar = ones(size_grid,1);
stateTransMat = ones(size_grid,size_grid);

% % Initialize covaraince and mean function
% covUpdateF_prev = basisCov;
% meanFuncF_prev = zeros(size_grid,1);
% meanFuncBasis = zeros(size_grid,1);

discTimeFiltCnst = sample_time/(tauFiltColor+sample_time);

%-------------------------------------------------------------------------
    desParaCur = [betaBarPrev AzLimPrev];
    responseCurr = perfmIndexVal;
    
    covCurPara = zeros(1,size_grid);
    for ii = 1:size_grid
        covCurPara(1,ii) = covFuncEval(desParaCur,basisParaGrid(ii,:),hypParaMat,sigmaSqKern);
    end 
    covCurrPoint = covFuncEval(desParaCur,desParaCur,hypParaMat,sigmaSqKern);
    
    covUpdateVecColored = covCurPara*invBasisCov;
        % Predicition step of Kalman filter
        meanFuncPredCurPt = meanFuncIC + covUpdateVecColored*(meanFuncF_prev-meanFuncIC*ones(size_grid,1));
        obsMat = covUpdateVecColored;

        I_mat = eye(size_grid);  
        stateTransMat = eye(size_grid);

        colNoisObsMat = covUpdateVecColored;
           
        colNoisTransMat = (1-discTimeFiltCnst);
        % Colored noise observation matrix
        Hprime = obsMat*stateTransMat - colNoisTransMat*colNoisObsMat;
        % Covariance at the current operating point
        covUpdatedP = covCurrPoint + covUpdateVecColored*(covUpdateF_prev - basisCov)*covUpdateVecColored';
        % Kalman Gain Calculation 
        kalGainMatrixNew = covUpdateF_prev*Hprime'*((Hprime*covUpdateF_prev*Hprime' + sigmaSq)^(-1));
        % Update step of Kalman filter (update states and covaraiance
        meanFuncFupdate = meanFuncF_prev + kalGainMatrixNew*(responseCurr - meanFuncPredCurPt);
%         covUpdatedF = (I_mat - kalGainMatrixNew*Hprime)*covUpdateF_prev*(I_mat - kalGainMatrixNew*Hprime)' + kalGainMatrixNew*sigmaSq*kalGainMatrixNew';
        covUpdatedF = (I_mat - kalGainMatrixNew*Hprime)*covUpdateF_prev*(I_mat - kalGainMatrixNew*Hprime)';

        
        predVar = diag(covUpdatedF);
        [J_hatStar,maxPerfmVal_index] = max(meanFuncFupdate);
        optBasisVec = basisParaGrid(maxPerfmVal_index,:);        
end 

%--------------------------------------------------------------------------    
function [covFuncValue] = covFuncEval(designPt1,designPt2,hypParaMat,alphaSq)
%     covFuncValue = noiseVar*exp(-1/(2*theta^2)*(abs(designPt1-designPt2))^2);
    covFuncValue = alphaSq*exp(-0.5*(designPt1-designPt2)*inv(hypParaMat)*(designPt1-designPt2)');
end 
    

