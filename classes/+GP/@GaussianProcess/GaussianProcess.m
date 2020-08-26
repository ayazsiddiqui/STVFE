classdef GaussianProcess
    %GAUSSIANPROCESS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        numberOfSpatialInputs
    end
    
    properties
       spatialKernel; 
    end
    
    %% constructor
    methods
        function obj = GaussianProcess(arg1)
            %GAUSSIANPROCESS Construct an instance of this class
            %   Detailed explanation goes here
            
            
            obj.spatialKernel = SquaredExponentialKernel;
            
            obj.numberOfSpatialInputs = arg1;
            
        end
        
    end
    %% setters
    methods
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    %% getters
    methods
        
    end
    
    %% other methods
    methods
        
    end
end

