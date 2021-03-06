function [hold,basVecMaxEntropyRGPC,infMetricStore,... 
    redBasisParaGridLogic] = ReduceColored_RGP_Ayaz(betaBarPrev,AzLimPrev,predVar,basisParaGrid,... 
    size_grid,zScoreRej,J_hat,redBasisParaGridLogicPrev)
    
    prevDesParam = [betaBarPrev AzLimPrev];
    
    redBasisParaGridLogic = redBasisParaGridLogicPrev;
    [~,maxJhatIndex] = max(J_hat);
    
    indexBasisVec = [1:1:size_grid]'; % This used to replace the find function.   
         
% Rejection here is done based error bars placed on mean function
% approximation at basis vector locations. The error bars are generated by:
% J_hat +/- 1.96*sigma which corresponds to 95%)

    meanFuncUpperErrorBar = J_hat + zScoreRej*sqrt(predVar);
    meanFuncLowerErrorBar = J_hat - zScoreRej*sqrt(predVar);
    
    meanFuncStarLowerErrorBar = meanFuncLowerErrorBar(maxJhatIndex);
    for ii = 1:size_grid
        if meanFuncUpperErrorBar(ii) >= meanFuncStarLowerErrorBar && redBasisParaGridLogicPrev(ii) == 1 
            redBasisParaGridLogic(ii) = 1;
        else 
            redBasisParaGridLogic(ii) = 0;
        end
    end
    
    % This is a patch to avoid taking the max or min of something with zero
    % length
    if sum(redBasisParaGridLogic(:,1)) <= 10
        redBasisParaGridLogic = redBasisParaGridLogicPrev;
    end 

    basisIndRedSpace = indexBasisVec.*redBasisParaGridLogic;
    minBasisVecInt = min(basisIndRedSpace);
        
    dist_update = zeros(size_grid,1);
    for c4 = 1:size_grid
        if redBasisParaGridLogic(c4,1) == 1
            dist_update(c4,1) = sqrt((prevDesParam(1,1)-basisParaGrid(c4,1))^2 + (prevDesParam(1,2)-basisParaGrid(c4,2))^2);
        else 
            dist_update(c4,1) = 0;
        end
    end 
        
    if sum(redBasisParaGridLogic(:,1)) <= 10 
        if minBasisVecInt < size_grid
            basVecMaxEntropyRGPC = basisParaGrid(minBasisVecInt+1,:);
        else 
            basVecMaxEntropyRGPC = basisParaGrid(minBasisVecInt,:);
        end 
        infMetricStore = zeros(size_grid,1);
    else 
        predVarReducedDesSpace = predVar.*redBasisParaGridLogic(:,1);
        [valMaxEntropy,basVecMaxEntropyIndex] = max(predVarReducedDesSpace.*exp(-75*dist_update.^2));

        basVecMaxEntropyRGPC = basisParaGrid(basVecMaxEntropyIndex,:);
        infMetricStore = predVar;
    end 

%--------------------------------------------------------------------------
% Convergence test
meanFuncConvTest = NaN(size_grid,1);
for mm = 1:size_grid
    if basisIndRedSpace(mm) > 0
        meanFuncConvTest(mm,1) = J_hat(mm,1);
    else
        meanFuncConvTest(mm,1) = NaN;
    end
end

maxMeanFuncEachStep = max(meanFuncConvTest);
minMeanFuncEachStep = min(meanFuncConvTest);
percentDiff = (maxMeanFuncEachStep - minMeanFuncEachStep)/(0.5*(maxMeanFuncEachStep+minMeanFuncEachStep));

if percentDiff <= 0.1 && sum(redBasisParaGridLogic(:,1)) < size_grid-1 
    hold = 1;
    [~,minConvIdx] = min(meanFuncConvTest);
    basVecMaxEntropyRGPC = basisParaGrid(minConvIdx,:);
else
    hold = 0;
end 