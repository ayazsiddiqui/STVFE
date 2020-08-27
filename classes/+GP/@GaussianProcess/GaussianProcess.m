classdef GaussianProcess
    %GAUSSIANPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        numberOfSpatialInputs
        spatialKernel
        temporalKernel
    end
    
    % hyper-parameter properties
    properties
        spatialLengthScale
        spatialCovAmp
        temporalLengthScale
        temporalCovAmp
        noiseVariance
    end
    
    properties (SetAccess = protected)
        kernelChoices = {'exponential','squaredExponential'};
    end
    
    %% constructor
    methods
        function obj = GaussianProcess(noSpatialIps,spatialKernel,...
                temporalKernel)
            % set number of spatial inputs
            obj.numberOfSpatialInputs = noSpatialIps;
            obj.temporalKernel        = 'exponential';
            obj.spatialKernel         = 'squaredExponential';
            
            if nargin>1
                % check for valid inputs
                obj.temporalKernel = temporalKernel;
                obj.spatialKernel  = spatialKernel;
            end
            
        end
        
    end
    %% setters
    methods
        % set spatial kernel
        function obj = set.spatialKernel(obj,value)
            obj.spatialKernel = obj.checkKernelValidity(value);
        end
        
        % set temporal kernel
        function obj = set.temporalKernel(obj,value)
            obj.temporalKernel = obj.checkKernelValidity(value);
        end
        
        % set spatial length scale
        function obj = set.spatialLengthScale(obj,value)
            obj.spatialLengthScale = value(:);
        end
        
    end
    
    %% getters
    methods
        
    end
    
    %% private methods. input checks and stuff
    methods (Access = protected)
        % check kernel choice validity
        function val = checkKernelValidity(obj,ipKernel)
            if ismember(ipKernel,obj.kernelChoices)
                val = ipKernel;
            else
                error(['Only ',repmat('%s, ',1,numel(obj.kernelChoices)-1),...
                    'and %s are valid entries for kernels.',...
                    ' You entered %s.'],obj.kernelChoices{:},ipKernel);
            end
            % set values
            switch ipKernel
                case 'exponential'
                    val = @ExponentialKernel;
                case 'squaredExponential'
                    val = @SquaredExponentialKernel;
            end
        end
        
    end
    
    %%  methods related to covariance calculations
    methods
        % calculate spatial covariance
        function val = calcSpatialCovariance(obj,x1,x2)
            val = obj.spatialKernel(x1,x2,[obj.spatialCovAmp;...
                obj.spatialLengthScale]);
        end
        
        % calculate temporal covariance
        function val = calcTemporalCovariance(obj,t1,t2)
            val = obj.temporalKernel(t1,t2,[obj.temporalCovAmp;...
                obj.temporalLengthScale]);
        end
        
        % calculate total covariance
        function val = calcTotalCovariance(obj,x1,x2,t1,t2)
            val = obj.calcSpatialCovariance(x1,x2)*...
                obj.calcTemporalCovariance(t1,t2);
        end
        
        % make spatial covariance matrix
        function val = makeSpatialCovarianceMatrix(obj,X)
            % local variables
            nPOints = size(X,2);
            % preallocate
            val = zeros(nPOints);
            % loop throught each point and make lower triu matrix
            for ii = 1:nPOints
                for jj = 1:ii
                    val(ii,jj) = obj.calcSpatialCovariance(X(:,ii),X(:,jj));
                end
            end
            % make the total matrix
            val = val + tril(val,-1)';
        end
        
        % make temporal covariance matrix
        function val = makeTemporalCovarianceMatrix(obj,T)
            % local variables
            nPoints = size(T,2);
            % preallocate
            val = zeros(nPoints);
            % loop throught each point and make lower triu matrix
            for ii = 1:nPoints
                for jj = 1:ii
                    val(ii,jj) = obj.calcTemporalCovariance(T(:,ii),T(:,jj));
                end
            end
            % make the total matrix
            val = val + tril(val,-1)';
        end
        
        % make total covariance matrix
        function val = makeTotalCovarianceMatrix(obj,XT)
            % local variables
            nPoints = size(XT,2);
            % preallocate
            val = zeros(nPoints);
            % loop throught each point and make lower triu matrix
            for ii = 1:nPoints
                for jj = 1:ii
                    val(ii,jj) = obj.calcTotalCovariance(XT(1:end-1,ii),...
                        XT(1:end-1,jj),XT(end,ii),XT(end,jj));
                end
            end
            % make the total matrix
            val = val + tril(val,-1)';
        end
    end
    
    %% methods to caclualate various marginal likelihoods
    methods (Access = protected)
        % calculate spatial marginal likelihood.
        function val = calcSpatialMarginalLikelihood(obj,X,y,hyp)
            % update the object hyper-parameters
            obj.spatialCovAmp = hyp(1);
            obj.spatialLengthScale = hyp(2:end);
            % local variables
            noiseVar = obj.noiseVariance;
            nPoints  = numel(y);
            Kf = obj.makeSpatialCovarianceMatrix(X);
            Ky = Kf + eye(nPoints)*noiseVar;
            % log likelihood equation. Source: Rasmussen pg. 113
            dataFitPenalty      = -0.5*y'*(Ky\y);
            complexityPenalty   = -0.5*log(det(Ky));
            normalizingConstant = -0.5*nPoints*log(2*pi);
            val = dataFitPenalty + complexityPenalty + normalizingConstant;
        end
        
        % calculate temporal marginal likelihood
        function val = calcTemporalMarginalLikelihood(obj,T,y,hyp)
            % update the object hyper-parameters
            obj.temporalCovAmp = hyp(1);
            obj.temporalLengthScale = hyp(2);
            % local variables
            noiseVar = obj.noiseVariance;
            nPoints  = numel(y);
            Kf = obj.makeTemporalCovarianceMatrix(T);
            Ky = Kf + eye(nPoints)*noiseVar;
            % log likelihood equation. Source: Rasmussen pg. 113
            dataFitPenalty      = -0.5*y'*(Ky\y);
            complexityPenalty   = -0.5*log(det(Ky));
            normalizingConstant = -0.5*nPoints*log(2*pi);
            val = dataFitPenalty + complexityPenalty + normalizingConstant;
        end
        
        % calculate total marginal likelihood
        function val = calcTotalMarginalLikelihood(obj,XT,y,hyp)
            % update the object hyper-parameters
            obj.spatialCovAmp       = hyp(1);
            obj.spatialLengthScale  = hyp(2:end-2);
            obj.temporalCovAmp      = hyp(end-1);
            obj.temporalLengthScale = hyp(end);
            % local variables
            noiseVar = obj.noiseVariance;
            nPoints  = numel(y);
            Kf = obj.makeTotalCovarianceMatrix(XT);
            Ky = Kf + eye(nPoints)*noiseVar;
            % log likelihood equation. Source: Rasmussen pg. 113
            dataFitPenalty      = -0.5*y'*(Ky\y);
            complexityPenalty   = -0.5*log(det(Ky));
            normalizingConstant = -0.5*nPoints*log(2*pi);
            val = dataFitPenalty + complexityPenalty + normalizingConstant;
        end
        
    end
    
    %% methods to optimize hyper-parameters
    methods
        % optimize spatial hyper-parameters
        function val = findOptSpatialHyperParams(obj,X,y,iniGuess)
            % set optimization solver options
            options = optimoptions('fmincon','Display','notify');
            % local variables
            nHyps = numel(iniGuess);
            lb = 0.001*ones(nHyps,1);
            % find vals
            optVals = ...
                fmincon(@(hyp)-obj.calcSpatialMarginalLikelihood(X,y,hyp),...
                iniGuess,[],[],[],[],lb,[],[],options);
            % output
            val.opt_spatialCovAmp = optVals(1);
            val.opt_spatialLengthScale = optVals(2:end);
        end
        
        % optimize temporal hyper-parameters
        function val = findOptTemporalHyperParams(obj,T,y,iniGuess)
            % set optimization solver options
            options = optimoptions('fmincon','Display','notify');
            % local variables
            nHyps = numel(iniGuess);
            lb = 0.001*ones(nHyps,1);
            % find vals
            optVals = ...
                fmincon(@(hyp)-obj.calcTemporalMarginalLikelihood(T,y,hyp),...
                iniGuess,[],[],[],[],lb,[],[],options);
            % output
            val.opt_temporalCovAmp = optVals(1);
            val.opt_temporalLengthScale = optVals(2);
        end
        
        % optimize spatio-temporal hyper-parameters
        function val = findOptSpatioTemporalHyperParams(obj,XT,y,iniGuess)
            % set optimization solver options
            options = optimoptions('fmincon','Display','notify');
            % local variables
            nHyps = numel(iniGuess);
            lb = 0.001*ones(nHyps,1);
            % find vals
            optVals = ...
                fmincon(@(hyp)-obj.calcTotalMarginalLikelihood(XT,y,hyp),...
                iniGuess,[],[],[],[],lb,[],[],options);
            % output
            val.opt_spatialCovAmp       = optVals(1);
            val.opt_spatialLengthScale  = optVals(2:end-2);
            val.opt_temporalCovAmp      = optVals(end-1);
            val.opt_temporalLengthScale = optVals(end);
        end
        
    end
    
    %%
    
    
    
    
    
end



