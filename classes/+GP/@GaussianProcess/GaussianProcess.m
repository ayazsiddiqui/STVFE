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
    methods (Access = protected,Hidden)
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
    
    %% other methods
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
        % make covariance matrix
        function val = makeSpatialCovarianceMatrix(obj,X1,X2)
            % local variables
            nPoints = size(X1,2);
            % preallocate
            val = zeros(nPoints);
            % loop throught each point and make lower triu matrix
            for ii = 1:nPoints
                for jj = 1:ii
                    val(ii,jj) = obj.calcSpatialCovariance(X1(:,ii),X2(:,ii));
                end
            end
            % make the total matrix
            val = val + triu(val,1)';
        end
            
        
    end
end

