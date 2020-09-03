classdef RecursiveGaussianProcess < GP.GaussianProcess
    
    % constructor properties
    properties (SetAccess = immutable)
        xBasis
    end
    
    %
    properties
        spatialCovMat
        meanFnVector
    end
    
    %% constructor
    methods
        function obj = RecursiveGaussianProcess(spatialKernel,...
                meanFn,xBasis)
            % call superclass constructor
            obj@GP.GaussianProcess(spatialKernel,'alwaysOne',meanFn);
            % set class properties
            obj.xBasis     = xBasis;
        end
        
    end
    
    %% initialization related methods
    methods
        
    end
    
    %% regression related methods
    methods
        
        function [predMean,postVarMat] = ...
                calcPredMeanAndPostVar(obj,muGt_1,cGt_1,xt,yt)            
            % extract values from basic vector
            xB = obj.xBasis;
            % number of design/training points
            noTP = size(xB,2);
            % calculate mean and covariance at candidate point
            mXt = obj.meanFunction(xt);
            kXtXt = obj.calcSpatialCovariance(xt,xt);
            % calculate covariance of candidate wrt design points
            kXtX = NaN(1,noTP);
            for ii = 1:noTP
                kXtX(1,ii) = obj.calcSpatialCovariance(xt,xB(:,ii));
            end
            % calculate Jt as per Huber Eqn. (8)
            Jt = kXtX/obj.spatialCovMat;
            % calculate B as per Huber Eqn. (7)
            B = kXtXt - Jt*kXtX';
            % calculate muP as per Huber Eqn. (6)
            muP = mXt + Jt*(muGt_1 - obj.meanFnVector)';
            % calculate cP as per Huber Eqn. (9)
            cP = B + Jt*cGt_1*Jt';
            % calculate Gt (kalman gain matrix) as per Huber Eqn. (12)
            Gt = cGt_1*Jt'*(cP + obj.noiseVariance)^-1;
            % calculate mean at step t as per Huber Eqn. (10)
            predMean = muGt_1' + Gt*(yt - muP);
            % calculate covariance at step t as per Huber Eqn. (11)
            postVarMat = cGt_1 - Gt*Jt*cGt_1;         
        end
        
    end
    
    
    
    
    
end

