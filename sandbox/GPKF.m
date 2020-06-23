classdef GPKF
    %GPKF classdef containing various functions for performing operations
    %using the Gassian process kalman filter
    %   gpkf = GPKF(n) where n is the number of spatial inputs
    
    %% properties
    properties
        noSpatialIps
        temporalKernel
        acquisitionFunction
    end
    
    %% methods
    methods
        
        % % % %         constructor
        function obj = GPKF(noSpatialIps,temporalKernel,acquisitionFunction)
            obj.noSpatialIps = noSpatialIps;
            % % % check validity of temporalKernel
            temporalKernelChoice = {'exponential','squaredExponential'};
            if ismember(temporalKernel,temporalKernelChoice)
                obj.temporalKernel = temporalKernel;
            else
                error(['Only ',repmat('%s, ',1,numel(temporalKernelChoice)-1),...
                    'and %s are valid entries for temporalKernel.',...
                    ' You entered %s.'],temporalKernelChoice{:},temporalKernel);
            end
            % % % check validity of acquisition function
            aquiFunChoice = {'upperConfidenceBound','expectedImprovement'};
            if ismember(acquisitionFunction,aquiFunChoice)
                obj.acquisitionFunction = acquisitionFunction;
            else
                error(['Only ',repmat('%s, ',1,numel(aquiFunChoice)-1),...
                    'and %s are valid entries for temporalKernel.',...
                    ' You entered %s.'],aquiFunChoice{:},acquisitionFunction);
            end
            
        end
        
        % % % %         mean function
        function val = meanFunction(obj,x)
            % % zero mean function
            val = 0*(x'*x)*obj.noSpatialIps;
        end
        
        % % % %         spatial kernel: squared exponential
        function val = spatialCovariance(obj,s1,s2,covAmp,lengthScales)
            % % covariance equations
            val = covAmp*exp(-0.5*(s1-s2)'*(eye(obj.noSpatialIps)...
                ./lengthScales.^2)*(s1-s2));
        end
        
        % % % %         form covariance matrix and mean vector
        function covMat = buildSpatialCovMat(obj,x,covAmp,lengthScales)
            % number of points = number of columns
            noPts = size(x,2);
            % preallocate matricx
            covMat = zeros(noPts);
            % form lower triangular covariance matrix for covMat
            for ii = 1:noPts
                for jj = ii:noPts
                    covMat(ii,jj) = obj.spatialCovariance(x(:,ii),x(:,jj),...
                        covAmp,lengthScales);
                end
            end
            % form the total covariance matrix
            covMat = covMat + triu(covMat,1)';
        end
        
        % % % %         temporal kernel
        function val = temporalCovariance(obj,t1,t2,timeScale)
            % % covariance equation
            switch obj.temporalKernel
                case 'exponential'
                    val = 1*exp(-abs(t2-t1)/timeScale);
                case 'squaredExponential'
                    val = 1*exp(-0.5*((t2-t1)^2)/timeScale^2);
            end
        end
        
        % % % %         acquisition function
        function val = calcAcquisitionFun(obj,predMean,postVar,fBest,varargin)
            % % parse inputs
            pp = inputParser;
            addParameter(pp,'explorationConstant',2,@(x) isnumeric(x));
            addParameter(pp,'exploitationConstant',1,@(x) isnumeric(x));
            parse(pp,varargin{:});
            
            switch obj.acquisitionFunction
                case 'upperConfidenceBound'
                    % calculate UCB
                    val = pp.Results.exploitationConstant*predMean + ...
                        2^(pp.Results.explorationConstant).*postVar;
                case 'expectedImprovement'
                    % calculate standard deviation
                    stdDev = postVar.^0.5;
                    % make standard normal distribution
                    stdNormDis = gmdistribution(0,1);
                    % calculate z
                    Z(stdDev>0,1) = (predMean(:) - fBest)./stdDev(:);
                    Z(stdDev<=0,1) = 0;
                    % calculate expected improvement
                    val = (predMean - fBest).*cdf(stdNormDis,Z) + ...
                        stdDev.*pdf(stdNormDis,Z);
                    val(stdDev<=0) = 0;
            end
            
        end
        
        % % % %         calculate covariance as a product of the two covariances
        function val = calcTotCovariance(obj,x1,x2,hyperParams)
            % % covariance amplitude or variance of latent function
            covAmp = hyperParams(1);
            % % length scales for spatial covariance
            lenScale = hyperParams(2:obj.noSpatialIps+1);
            % % time scale
            timeScale = hyperParams(obj.noSpatialIps+2);
            % % k = k_s*k_t
            val = obj.spatialCovariance(x1(1:end-1),x2(1:end-1),covAmp,lenScale)...
                *obj.temporalCovariance(x1(end),x2(end),timeScale);
        end
        
        % % % %         form covariance matrix and mean vector
        function [covMat,meanVec] = buildCovMatAndMeanVec(obj,x,hyperParams)
            % number of points = number of columns
            noPts = size(x,2);
            % initial matrices
            covMat = zeros(noPts);
            meanVec = NaN(noPts,1);
            % form lower triangular covariance matrix for covMat and
            % calculate mean
            for ii = 1:noPts
                for jj = ii:noPts
                    covMat(ii,jj) = obj.calcTotCovariance(x(:,ii),x(:,jj),...
                        hyperParams);
                end
                meanVec(ii,1) = obj.meanFunction(x(:,ii));
            end
            % form the total covariance matrix
            covMat = covMat + triu(covMat,1)';
        end
        
        % % % %         calculate marginal likelihood
        function logP = calcMarginalLikelihood(obj,x,y,hyperParams)
            % determine number of training points
            noTP = size(x,2);
            % build the covariance matrix
            kX = obj.buildCovMatAndMeanVec(x,hyperParams);
            % add signal noise to the covariance matrix
            kX = kX + eye(noTP)*hyperParams(obj.noSpatialIps + 3);
            % the data fit part
            dataFit = 0.5*y'*(kX\y);
            % complexity penalty
            complexPen = 0.5*log(det(kX));
            % normalizing constant
            normConstant = 0.5*noTP*log(2*pi);
            % marginal likelihood
            logP = - dataFit - complexPen - normConstant;
        end
        
        % % % %         GPKF initialization
        function val = GpkfInitialize(obj,xDomain,timeScale,...
                timeStep,varargin)
            % % parse inputs
            pp = inputParser;
            addParameter(pp,'approximationOrder',6,...
                @(x)assert(isnumeric(x) && ~mod(x,2),...
                'Value must be even.'));
            parse(pp,varargin{:});
            N = pp.Results.approximationOrder;
            
            % % total number of points in the entire domain of interest
            xDomainNP = size(xDomain,2);
            
            switch obj.temporalKernel
                case 'exponential'
                    % % calculate F,H,Q as per Carron Eqn. (14)
                    F = exp(-timeStep/timeScale);
                    H = sqrt(2/timeScale);
                    G = 1;
                    Q = (1 - exp(-2*timeStep/timeScale))/(2/timeScale);
                    % % solve the Lyapunov equation for X
                    sigma0 = lyap(F,G*G');
                    % % outputs
                    val.Amat = eye(xDomainNP)*F;
                    val.Hmat = eye(xDomainNP)*H;
                    val.Qmat = eye(xDomainNP)*Q;
                    val.sig0Mat = eye(xDomainNP)*sigma0;
                    val.s0 = zeros(xDomainNP,1);
                    
                case 'squaredExponential'
                    % % find the transfer function as per the Hartinkainen paper
                    syms x
                    px = 0;
                    % % Hartinkainen paper Eqn. (11)
                    for n = 0:N
                        px = px + ((x^(2*n))*factorial(N)*((-1)^n)*...
                            (2/(timeScale^2))^(N-n))...
                            /factorial(n);
                    end
                    % % find the roots of the above polynomial
                    rts = vpasolve(px,x);
                    % % locate the roots with negative real parts
                    negReal = rts(real(rts) < 0);
                    % % make transfer function out of the negative real parts roots
                    H_iw = vpa(expand(prod(x-negReal)));
                    % % find the coefficients of the polynomial
                    coEffs = coeffs(H_iw,x);
                    % break the coefficients in real and imaginary parts and
                    % eliminate numbers lower than eps
                    nRound = 8;
                    coEffs = obj.removeEPS(coEffs,nRound);
                    % % normalize them by dividing by the highest degree
                    % % coefficient
                    coEffs = coEffs./coEffs(end);
                    % % form the F, G, and H matrices as per Carron Eqn. (8)
                    F = [zeros(N-1,1) eye(N-1); -coEffs(1:end-1)];
                    G = [zeros(N-1,1);1];
                    % % calculate the numerator
                    b0 = sqrt(factorial(N)*((2/(timeScale^2))^N)...
                        *sqrt(pi*2*timeScale^2));
                    H = [b0 zeros(1,N-1)];
                    sigma0 = obj.removeEPS(lyap(F,G*G'),nRound);
                    % % calculate the discretized values
                    syms tau
                    % % use cayley hamilton theorem to calcualte e^Ft
                    eFt = obj.cayleyHamilton(F);
                    % % calculate Fbar using the above expression
                    Fbar = obj.removeEPS(subs(eFt,tau,timeStep),nRound);
                    % % evaluate Qbar, very computationally expensive
                    Qsym = eFt*(G*G')*eFt';
                    Qint = NaN(N);
                    for ii = 1:N^2
                        fun = matlabFunction(Qsym(ii));
                        Qint(ii) = integral(fun,0,timeStep);
                    end
                    % % remove numbers lower than eps
                    Qbar = obj.removeEPS(Qint,nRound);
                    % % outputs
                    % initialize matrices as cell matrices
                    Amat = cell(xDomainNP);
                    Hmat = cell(xDomainNP);
                    Qmat = cell(xDomainNP);
                    sig0Mat = cell(xDomainNP);
                    
                    % form the block diagonal matrices
                    for ii = 1:xDomainNP
                        for jj = 1:xDomainNP
                            if ii == jj
                                Amat{ii,jj} = Fbar;
                                Hmat{ii,jj} = H;
                                Qmat{ii,jj} = Qbar;
                                sig0Mat{ii,jj} = sigma0;
                            else
                                Amat{ii,jj} = zeros(N);
                                Hmat{ii,jj} = zeros(1,N);
                                Qmat{ii,jj} = zeros(N);
                                sig0Mat{ii,jj} = zeros(N);
                            end
                        end
                    end
                    % convert them to matrices and send to output structure
                    val.Amat = cell2mat(Amat);
                    val.Hmat = cell2mat(Hmat);
                    val.Qmat = cell2mat(Qmat);
                    val.sig0Mat = cell2mat(sig0Mat);
                    val.s0 = zeros(xDomainNP*N,1);
            end
        end
        
        % % % %         Cayley Hamilton theorem implementation
        function eAt = cayleyHamilton(obj,A)
            % order of matrix
            n = length(A);
            % eigen values
            eVals = eig(A);
            reVals = real(eVals);
            ieVals = imag(eVals);
            % define t
            syms tau
            lhs = sym('lhs',[n,1]);
            rhs = NaN(n,n);
            % populate the LHS and RHS matrices
            for ii = 1:n
                lhs(ii) = exp(reVals(ii)*tau)*(cos(abs(ieVals(ii))*tau) + ...
                    1i*sign(ieVals(ii))*sin(abs(ieVals(ii))*tau));
                for jj = 1:n
                    rhs(ii,jj) = eVals(ii)^(jj-1);
                end
            end
            % solve for alpha
            alp = simplify(inv(rhs)*lhs(:));
            eAt = zeros(n);
            % calculate the e^At matrix
            for ii = 1:n
                eAt = alp(ii).*A^(ii-1) + eAt;
            end
            % simplify the symbolic expression
            eAt = vpa(simplify(eAt),4);
        end
        
        % % % %         remove values lower than eps from arrays
        function val = removeEPS(obj,xx,nRound)
            % eliminate numbers lower than eps
            realParts = double(real(xx));
            realParts(realParts <= eps) = 0;
            imagParts = double(imag(xx));
            imagParts(imagParts <= eps) = 0;
            % round up to nRound decimal places
            val = round(realParts,nRound) + 1i*round(imagParts,nRound);
        end
        
        % % % %         Kalman estimation as per jp Algorithm 1
        function [F_t,sigF_t,skp1_kp1,ckp1_kp1,varargout] = ...
                gpkfKalmanEstimation(obj,xMeasure,sk_k,ck_k,Mk,yk,...
                Ks_12,Amat,Qmat,Hmat,noiseVar)
            % % number of measurable points which is subset of xDomain
            xMeasureNP = size(xMeasure,2);
            % % number of points visited at each step which is a subset of xMeasure
            MkNP = size(Mk,2);
            % % R matrix as per Carron conf. paper Eqn. (12)
            Rmat = eye(MkNP)*noiseVar;
            % % indicator matrix to find which points are visited at each iteration
            Ik = zeros(MkNP,xMeasureNP);
            % % find points from xMeasure visited at iteration k
            lia = ismember(xMeasure',Mk','rows');
            lia = find(lia); % covert to numerical array instead of logical
            % % populate the Ik matrix
            for ii = 1:MkNP
                Ik(ii,lia(ii)) = 1;
            end
            % % C matrix as per Carron conf. paper Eqn. (12)
            Cmat = Ik*Ks_12*Hmat;
            % % Kalman filter equations as per Carron conf. paper Eqn. (6)
            skp1_k = Amat*sk_k; % Eqn. (6a)
            ckp1_k = Amat*ck_k*Amat' + Qmat; % Eqn. (6b)
            Lkp1 = ckp1_k*Cmat'/(Cmat*ckp1_k*Cmat' + Rmat); % Eqn (6e)
            % set kalman gain to zero if yk is empty
            if isempty(yk)
                skp1_kp1 = skp1_k; % Eqn (6c)
            else
                skp1_kp1 = skp1_k + Lkp1*(yk - Cmat*skp1_k); % Eqn (6c)
            end
            ckp1_kp1 = ckp1_k - Lkp1*Cmat*ckp1_k; % Eqn (6d)
            % % process estimate and covariance as per Todescato algortihm 1
            F_t = Ks_12*Hmat*skp1_kp1; % Eqn. (13)
            sigF_t = Ks_12*Hmat*ckp1_kp1*Hmat'*Ks_12;
            % % varibale outputs
            varargout{1} = Ik;
        end
        
        % % % %         Regression as per jp section 5
        function [predMean,postVar] = gpkfRegression(obj,xDomain,xPredict,...
                F_t,sigF_t,Ks,hyperParam)
            % % number of points in the discretized domain
            xDomainNP = size(xDomain,2);
            % % number of points over which we want to acquire predictions
            xPredictNp = size(xPredict,2);
            % % extract values from hyper parameters
            covAmp = hyperParam(1);
            lengthScales = hyperParam(2:end-2);
            timeScale = hyperParam(end-1);
            % % Regression as per section 5 of Todescato journal paper
            % % temporations kernel value at tau = 0
            h0 = obj.temporalCovariance(0,0,timeScale);
            % % multiply the spatial covariance matrix by h0
            Vf = h0*Ks;
            % % preallocate matrices
            sigmaX = NaN(xPredictNp,xDomainNP);
            Vx = NaN(xPredictNp,1);
            predMean = NaN(xPredictNp,1);
            postVar = NaN(xPredictNp,1);
            % % perform regression on each point in the domain
            for ii = 1:xPredictNp
                for jj = 1:xDomainNP
                    sigmaX(ii,jj) = h0*obj.spatialCovariance(xPredict(:,ii),xDomain(:,jj)...
                        ,covAmp,lengthScales);
                end
                Vx(ii,1) = h0*obj.spatialCovariance(xPredict(:,ii),xPredict(:,ii)...
                    ,covAmp,lengthScales);
                % % predicted mean as per Todescato Eqn. (17)
                predMean(ii,1) = sigmaX(ii,:)*(Vf\F_t);
                % % posterior variance as per Todescato Eqn. (18)
                postVar(ii,:) = Vx(ii,1) - sigmaX(ii,:)*(Vf\(Vf - sigF_t))*...
                    (Vf\sigmaX(ii,:)');
            end
            postVar = obj.removeEPS(postVar,6);
        end
        
        % % % %         future predictions using GPKF
        function op = predictionGPKF(obj,xMeasure,sk_k,ck_k,Mk,yk,...
                Ks_12,Amat,Qmat,Hmat,xPredict,Ks,hyperParam,predHorizon)
            % % % preallocate matrices
            F_t =  NaN(size(xMeasure,2),predHorizon);
            sigF_t = NaN(size(xMeasure,2),size(xMeasure,2),predHorizon);
            predMean = NaN(size(xPredict,2),predHorizon);
            postVar =  NaN(size(xPredict,2),predHorizon);
            stdDev =  NaN(size(xPredict,2),predHorizon);
            upperBound = NaN(size(xPredict,2),predHorizon);
            lowerBound = NaN(size(xPredict,2),predHorizon);
            predMeanAtLoc = NaN(1,predHorizon);
            postVarAtLoc = NaN(1,predHorizon);
            % % % obtain prediction mean and posterior variance over the
            % % % prediction horizon
            for jj = 1:predHorizon
                if jj == 1
                    ykPassed = yk;
                else
                    ykPassed = [];
                end
                
                % % % stepwise update of kalman state estimate and covariance
                [F_t(:,jj),sigF_t(:,:,jj),skp1_kp1,ckp1_kp1,Ik] = ...
                    obj.gpkfKalmanEstimation(xMeasure,sk_k,ck_k,Mk(jj),...
                    ykPassed,Ks_12,Amat,Qmat,Hmat,...
                    hyperParam(end));
                % % % kalman regression to calculate mean and variance
                [predMean(:,jj),postVar(:,jj)] = ...
                    obj.gpkfRegression(xMeasure,xPredict,...
                    F_t(:,jj),sigF_t(:,:,jj),Ks,hyperParam);
                % % % store mean and variance at current location for later
                predMeanAtLoc(jj) = predMean(logical(Ik'),jj);
                postVarAtLoc(jj) = postVar(logical(Ik'),jj);
                % % % remove real or imaginary parts lower than eps
                stdDev(:,jj) = sqrt(obj.removeEPS(postVar(:,jj),5));
                % % % upper bounds = mean + x*(standard deviation)
                upperBound(:,jj) = predMean(:,jj) + 1*stdDev(:,jj);
                % % % lower bounds = mean + x*(standard deviation)
                lowerBound(:,jj) = predMean(:,jj) - 1*stdDev(:,jj);
                % % % update previous step information
                sk_k = skp1_kp1;
                ck_k = ckp1_kp1;
            end
            
            op.predMeanPrediction = predMean;
            op.postVarPrediction = postVar;
            op.predMeanAtLoc = predMeanAtLoc;
            op.postVarAtLoc = postVarAtLoc;
        end
        
        % % % %         control using mpc
        function val = gpkfMPC(obj,xMeasure,sk_k,ck_k,Mk,yk,Ks_12,Amat,Qmat,...
                Hmat,xPredict,Ks,hyperParam,uAllowable,predHorizon,varargin)
            % % % parse input
            pp = inputParser;
            addParameter(pp,'explorationConstant',2,@(x) isnumeric(x));
            addParameter(pp,'exploitationConstant',1,@(x) isnumeric(x));
            parse(pp,varargin{:});
            % % % determine all possible control
            ctrlComb = makeBruteForceCombinations(uAllowable,predHorizon);
            % % % number of state trajectories
            nTraject = size(ctrlComb,1);
            % % % determine state trajectories, gpkf predictions,
            % % % and acquisition function values for all ctrlComb
            stateTrjectories = NaN(nTraject,predHorizon+1);
            stateTrjectories(:,1) = Mk;
            aqFunVal = NaN(nTraject,predHorizon+1);
            % % % start the looping
            for ii = 1:nTraject
                for jj = 2:predHorizon+1
                    % use previous state for following steps
                    stateTrjectories(ii,jj) = stateTrjectories(ii,jj-1)+...
                        ctrlComb(ii,jj-1);
                end
                % ensure the trajectory remains with bounds
                belowBounds = stateTrjectories(ii,:)<xMeasure(1);
                stateTrjectories(ii,belowBounds) = xMeasure(1);
                ctrlComb(ii,belowBounds(2:end)) = xMeasure(1) - ...
                    stateTrjectories(ii,belowBounds(2:end));
                
                aboveBounds = stateTrjectories(ii,:)>xMeasure(end);
                stateTrjectories(ii,aboveBounds) = xMeasure(end);
                ctrlComb(ii,aboveBounds(2:end)) = xMeasure(end) - ...
                    stateTrjectories(ii,aboveBounds(2:end));
            end
            % remove repeated state trajectories
            [stateTrjectories,ia,~] = unique(stateTrjectories,'rows');
            ctrlComb = ctrlComb(ia,:);
            for ii = 1:size(stateTrjectories,1)
                % obtain prediction mean and posterior variance over the
                % prediction horizon for each state trajectory
                gpkfPred = obj.predictionGPKF(xMeasure,sk_k,ck_k,...
                    stateTrjectories(ii,:),yk,Ks_12,Amat,Qmat,...
                    Hmat,xPredict,Ks,hyperParam,predHorizon+1);
                % % % calculate acquisition function values
                aqFunVal(ii,:) = obj.calcAcquisitionFun(gpkfPred.predMeanAtLoc(:),...
                    gpkfPred.postVarAtLoc(:),yk,'explorationConstant',...
                    pp.Results.explorationConstant,...
                    'exploitationConstant',pp.Results.exploitationConstant)';
            end
            % % % sum the acquisition function values along the column
            % direction
            sumAqFun = sum(aqFunVal,2);
            % find the max value and the corresponding state trajectory
            [maxVal,bestTrajIdx] = max(sumAqFun);
            % output these two
            val.optStateTrajectory = stateTrjectories(bestTrajIdx,:);
            val.optCtrlSeq = ctrlComb(bestTrajIdx,:);
            val.objFunVal = maxVal;
            
        end
        
        % % % %         traditional GP regression
        function [predMean,postVar] = traditionalGpRegression(obj,xVisited,yVisited,...
                xPredict,tPredict,hyperParams)
            % % number of points over which we want to acquire predictions
            xPredictNp = size(xPredict,2);
            % % total number of points visited
            xVisitedNp = size(xVisited,2);
            % % form the covariance matrix and mean vector
            [covMat,~] = obj.buildCovMatAndMeanVec(xVisited,hyperParams);
            % % add noise to the covariance
            covMat = covMat + eye(size(covMat))*hyperParams(end);
            % % preallocate matrices
            predMean = NaN(xPredictNp,1);
            postVar = NaN(xPredictNp,1);
            Kxstar_x = NaN(xPredictNp,xVisitedNp);
            % % prediction points
            xPredictNew = [xPredict;ones(1,xPredictNp)*tPredict];
            % % calculate prediction mean and covariance
            for ii = 1:xPredictNp
                for jj = 1:xVisitedNp
                    Kxstar_x(ii,jj) = obj.calcTotCovariance(xPredictNew(:,ii)...
                        ,xVisited(:,jj),hyperParams);
                end
                predMean(ii,1) = Kxstar_x(ii,:)*(covMat\yVisited);
                postVar(ii,1) = obj.calcTotCovariance(xPredictNew(:,ii),...
                    xPredictNew(:,ii),hyperParams) - Kxstar_x(ii,:)*...
                    (covMat\Kxstar_x(ii,:)');
            end
            
        end
        
    end
end

