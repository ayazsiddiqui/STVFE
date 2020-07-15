classdef GPKF
    %GPKF classdef containing various functions for performing operations
    %using the Gassian process kalman filter
    %   gpkf = GPKF(n) where n is the number of spatial inputs
    
    %% properties
    properties (SetAccess = private,Hidden)
        % constructore parameters
        p_noSpatialIps
        p_temporalKernel
        p_acquisitionFunction
    end
    
    properties
        % hyper-parameters
        p_spatialCovarianceAmp
        p_spatialLenghtScale
        p_temporalLenghtScale
        p_noiseVariance
        % gpkf parameters
        p_xMeasure
        p_gpfkTimeStep
        p_squaredExpApproxOrder = 2
        % acquisitioin function parameters
        p_explorationConstant
        p_exploitationConstant
    end
    
    properties (SetAccess = private)
        p_spatialCovMat
        p_spatialCovRoot
        p_gpfkInitVals
    end
    
    %% methods
    methods
        
        % % % %         constructor
        function obj = GPKF(noSpatialIps,temporalKernel,acquisitionFunction)
            obj.p_noSpatialIps = noSpatialIps;
            % % % check validity of temporalKernel
            temporalKernelChoice = {'exponential','squaredExponential'};
            if ismember(temporalKernel,temporalKernelChoice)
                obj.p_temporalKernel = temporalKernel;
            else
                error(['Only ',repmat('%s, ',1,numel(temporalKernelChoice)-1),...
                    'and %s are valid entries for temporalKernel.',...
                    ' You entered %s.'],temporalKernelChoice{:},temporalKernel);
            end
            % % % check validity of acquisition function
            aquiFunChoice = {'upperConfidenceBound','expectedImprovement'};
            if ismember(acquisitionFunction,aquiFunChoice)
                obj.p_acquisitionFunction = acquisitionFunction;
            else
                error(['Only ',repmat('%s, ',1,numel(aquiFunChoice)-1),...
                    'and %s are valid entries for temporalKernel.',...
                    ' You entered %s.'],aquiFunChoice{:},acquisitionFunction);
            end
            
        end
        
        % % % %         spatial covariance matrix
        function obj = m_CalcSpatialCovMat(obj)
            obj.p_spatialCovMat= obj.m_buildSpatialCovMat(obj.p_xMeasure);
        end
        
        % % % %         root of spatial covariance
        function obj = m_CalcSpatialCovRoot(obj)
            obj.p_spatialCovRoot =  sqrtm(obj.p_spatialCovMat);
        end
        
        % % % %         gpkf initialization values
        function obj =  gpfkInitialize(obj)
            obj.p_gpfkInitVals = obj.m_gpkfInitialize();
        end
        
        % % % %         mean function
        function val = m_meanFunction(obj,x)
            % % zero mean function
            val = 0*(x'*x)*obj.p_noSpatialIps;
        end
        
        % % % %         spatial kernel: squared exponential
        function val = m_CalcSpatialCovariance(obj,s1,s2)
            % % covariance equations
            val = obj.p_spatialCovarianceAmp*...
                exp(-0.5*(s1-s2)'*(eye(obj.p_noSpatialIps)...
                ./obj.p_spatialLenghtScale.^2)*(s1-s2));
        end
        
        % % % %         form covariance matrix
        function covMat = m_buildSpatialCovMat(obj,x)
            % number of points = number of columns
            noPts = size(x,2);
            % preallocate matricx
            covMat = zeros(noPts);
            % form lower triangular covariance matrix for covMat
            for ii = 1:noPts
                for jj = ii:noPts
                    covMat(ii,jj) = obj.m_CalcSpatialCovariance(x(:,ii),...
                        x(:,jj));
                end
            end
            % form the total covariance matrix
            covMat = covMat + triu(covMat,1)';
        end
        
        % % % %         temporal kernel
        function val = m_temporalCovariance(obj,t1,t2)
            % % covariance equation
            switch obj.p_temporalKernel
                case 'exponential'
                    val = 1*exp(-abs(t2-t1)/obj.p_temporalLenghtScale);
                case 'squaredExponential'
                    val = 1*exp(-0.5*((t2-t1)^2)/obj.p_temporalLenghtScale^2);
            end
        end
        
        % % % %         acquisition function
        function val = m_calcAcquisitionFun(obj,predMean,postVar,fBest)
            
            switch obj.p_acquisitionFunction
                case 'upperConfidenceBound'
                    % calculate UCB
                    val = obj.p_exploitationConstant.*predMean + ...
                        2^(obj.p_explorationConstant).*postVar;
                case 'expectedImprovement'
                    % calculate standard deviation
                    stdDev = postVar(:).^0.5;
                    predMean = predMean(:);
                    % make standard normal distribution
                    stdNormDis = gmdistribution(0,1);
                    % calculate z
                    Z(stdDev>0,1) = (predMean(stdDev>0) - fBest)...
                        ./stdDev(stdDev>0);
                    Z(stdDev<=0,1) = 0;
                    % calculate expected improvement
                    val = (predMean - fBest).*cdf(stdNormDis,Z) + ...
                        stdDev.*pdf(stdNormDis,Z);
                    val(stdDev<=0) = 0;
            end
        end
        
        % % % %         calculate covariance as a product of the two covariances
        function val = m_calcTotCovariance(obj,x1,x2)
            val = obj.m_CalcSpatialCovariance(x1(1:end-1),x2(1:end-1))...
                *obj.m_temporalCovariance(x1(end),x2(end));
        end
        
        % % % %         form covariance matrix and mean vector
        function [covMat,meanVec] = m_buildCovMatAndMeanVec(obj,x)
            % number of points = number of columns
            noPts = size(x,2);
            % initial matrices
            covMat = zeros(noPts);
            meanVec = NaN(noPts,1);
            % form lower triangular covariance matrix for covMat and
            % calculate mean
            for ii = 1:noPts
                for jj = ii:noPts
                    covMat(ii,jj) = obj.m_calcTotCovariance(x(:,ii),x(:,jj));
                end
                meanVec(ii,1) = obj.m_meanFunction(x(:,ii));
            end
            % form the total covariance matrix
            covMat = covMat + triu(covMat,1)';
        end
        
        % % % %         calculate marginal likelihood
        function logP = m_calcMarginalLikelihood(obj,x,y)
            % determine number of training points
            noTP = size(x,2);
            % build the covariance matrix
            kX = obj.m_buildCovMatAndMeanVec(x);
            % add signal noise to the covariance matrix
            kX = kX + eye(noTP)*obj.p_noiseVariance;
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
        function val = m_gpkfInitialize(obj)
            % % parse inputs
            N = obj.p_squaredExpApproxOrder;
            % % time step
            timeStep = obj.p_gpfkTimeStep;
            % % total number of points in the entire domain of interest
            xDomainNP = size(obj.p_xMeasure,2);
            % % switch cases
            switch obj.p_temporalKernel
                case 'exponential'
                    % % calculate F,H,Q as per Carron Eqn. (14)
                    F = exp(-timeStep/obj.p_temporalLenghtScale);
                    H = sqrt(2/obj.p_temporalLenghtScale);
                    G = 1;
                    Q = (1 - exp(-2*timeStep/obj.p_temporalLenghtScale))...
                        /(2/obj.p_temporalLenghtScale);
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
                            (2/(obj.p_temporalLenghtScale^2))^(N-n))...
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
                    b0 = sqrt(factorial(N)*((2/(obj.p_temporalLenghtScale^2))^N)...
                        *sqrt(pi*2*obj.p_temporalLenghtScale^2));
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
                m_gpkfKalmanEstimation(obj,sk_k,ck_k,Mk,yk)
            % % dummy variables
            xMeasure = obj.p_xMeasure;
            Ks_12 = obj.p_spatialCovRoot;
            Amat = obj.p_gpfkInitVals.Amat;
            Qmat = obj.p_gpfkInitVals.Qmat;
            Hmat = obj.p_gpfkInitVals.Hmat;
            % % number of measurable points which is subset of xDomain
            xMeasureNP = size(xMeasure,2);
            % % number of points visited at each step which is a subset of xMeasure
            MkNP = size(Mk,2);
            % % R matrix as per Carron conf. paper Eqn. (12)
            Rmat = eye(MkNP)*obj.p_noiseVariance;
            % % indicator matrix to find which points are visited at each iteration
            Ik = zeros(MkNP,xMeasureNP);
            % % populate the Ik matrix
            for ii = 1:MkNP
                distFromXm = Mk(ii) - xMeasure;
                distFromXm = round(distFromXm);
                [minDis,minIdx] = min(abs(distFromXm));
                if distFromXm(minIdx) > 0
                    Ik(ii,minIdx+1) = minDis/...
                        (xMeasure(minIdx+1)-xMeasure(minIdx));
                    Ik(ii,minIdx) = 1 - Ik(ii,minIdx+1);
                elseif distFromXm(minIdx) < 0
                    Ik(ii,minIdx-1) = minDis/...
                        (xMeasure(minIdx)-xMeasure(minIdx-1));
                    Ik(ii,minIdx) = 1 - Ik(ii,minIdx-1);
                else
                    Ik(ii,minIdx) = 1;
                end
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
        function [predMean,postVar] = m_gpkfRegression(obj,xPredict,...
                F_t,sigF_t)
            % % number of points in the discretized domain
            xDomainNP = size(obj.p_xMeasure,2);
            % % number of points over which we want to acquire predictions
            xPredictNp = size(xPredict,2);
            % % Regression as per section 5 of Todescato journal paper
            % % temporations kernel value at tau = 0
            h0 = obj.m_temporalCovariance(0,0);
            % % multiply the spatial covariance matrix by h0
            Vf = h0*obj.p_spatialCovMat;
            % % preallocate matrices
            sigmaX = NaN(xPredictNp,xDomainNP);
            Vx = NaN(xPredictNp,1);
            predMean = NaN(xPredictNp,1);
            postVar = NaN(xPredictNp,1);
            % % perform regression on each point in the domain
            for ii = 1:xPredictNp
                for jj = 1:xDomainNP
                    sigmaX(ii,jj) = h0*obj.m_CalcSpatialCovariance(...
                        xPredict(:,ii),obj.p_xMeasure(:,jj));
                end
                Vx(ii,1) = h0*obj.m_CalcSpatialCovariance(...
                    xPredict(:,ii),xPredict(:,ii));
                % % predicted mean as per Todescato Eqn. (17)
                predMean(ii,1) = sigmaX(ii,:)*(Vf\F_t);
                % % posterior variance as per Todescato Eqn. (18)
                postVar(ii,:) = Vx(ii,1) - sigmaX(ii,:)*(Vf\(Vf - sigF_t))*...
                    (Vf\sigmaX(ii,:)');
            end
            postVar = obj.removeEPS(postVar,6);
        end
        
        % % % %         future predictions using GPKF
        function op = m_predictionGPKF(obj,sk_k,ck_k,Mk,yk,traject,...
                predHorizon)
            
            % % % dummy variables
            matLengths = predHorizon + 1;
            % % % preallocate matrices
            predMeanAtLoc = NaN(1,matLengths);
            postVarAtLoc = NaN(1,matLengths);
            % % % obtain prediction mean and posterior variance over the
            % % % prediction horizon
            for jj = 1:matLengths
                if jj == 1
                    ykPassed = yk;
                    xLocation = Mk;
                else
                    ykPassed = [];
                    xLocation = traject(jj-1);
                end
                % % % stepwise update of kalman state estimate and covariance
                [F_t,sigF_t,skp1_kp1,ckp1_kp1] = ...
                    obj.m_gpkfKalmanEstimation(sk_k,ck_k,xLocation,ykPassed);
                % % % store mean and variance at current location for later
                [predMeanAtLoc(jj),postVarAtLoc(jj)]...
                    = obj.m_gpkfRegression(xLocation,F_t,sigF_t);
                % % % update previous step information
                sk_k = skp1_kp1;
                ck_k = ckp1_kp1;
            end
            
            op.predMeanAtLoc = predMeanAtLoc(2:end);
            op.postVarAtLoc = postVarAtLoc(2:end);
        end
        
        % % % %         control using mpc
        function val = m_gpkfMPC_bruteForce(obj,sk_k,ck_k,Mk,yk,...
                uAllowable,predHorizon)
            % % % dummy variables
            xMeasure = obj.p_xMeasure;
            % % % determine all possible control
            ctrlComb = makeBruteForceCombinations(uAllowable,predHorizon);
            % % % number of state trajectories
            nTraject = size(ctrlComb,1);
            % % % determine state trajectories, gpkf predictions,
            % % % and acquisition function values for all ctrlComb
            stateTrjectories = NaN(nTraject,predHorizon+1);
            stateTrjectories(:,1) = Mk;
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
            aqFunVal = NaN(numel(ia),predHorizon);

            for ii = 1:size(stateTrjectories,1)
                % obtain prediction mean and posterior variance over the
                % prediction horizon for each state trajectory
                gpkfPred = obj.m_predictionGPKF(sk_k,ck_k,...
                    Mk,yk,stateTrjectories(ii,2:end),...
                    predHorizon);
                % % % calculate acquisition function values
                aqFunVal(ii,:) = obj.m_calcAcquisitionFun(...
                    gpkfPred.predMeanAtLoc(:),...
                    gpkfPred.postVarAtLoc(:),yk)';
                % % %  weight the aqFun
                aqFunVal(ii,:) = (predHorizon:-1:1).*aqFunVal(ii,:);
            end
            % % % sum the acquisition function values along the column
            % direction
            sumAqFun = sum(aqFunVal,2);
            % find the max value and the corresponding state trajectory
            [maxVal,bestTrajIdx] = max(sumAqFun);
            % output these two
            val.optStateTrajectory = stateTrjectories(bestTrajIdx,2:end);
            val.optCtrlSeq = ctrlComb(bestTrajIdx,:);
            val.objFunVal = maxVal;
            
            fprintf('Exhaustive search optimization complete.\n')
        end
        
        function val = m_gpkfMPC_fmincon(obj,sk_k,ck_k,Mk,yk,...
                uLimits,predHorizon,numStarts)
            % % % set fmincon parameters
            % % % constraints
            % % % upper bounds
            ub = Mk + (1:predHorizon).*max(uLimits);
            ub(ub>=obj.p_xMeasure(end)) = obj.p_xMeasure(end);
            % % % lower bounds
            lb = Mk + (1:predHorizon).*min(uLimits);
            lb(lb<=obj.p_xMeasure(1)) = obj.p_xMeasure(1);
            % % % control bounds
            A3 = zeros(predHorizon-1,predHorizon);
            for ii = 1:predHorizon-1
                A3(ii,ii) = -1;
                A3(ii,ii+1) = 1;
            end
            b3 = ones(predHorizon-1,1)*max(uLimits);
            A4 = -A3;
            b4 = ones(predHorizon-1,1)*min(uLimits);
            A3 = [];
            A4 = [];
            b3 = [];
            b4 = [];
            % % % fmincon options
            options = optimoptions('fmincon','algorithm','sqp',...
                'Display','notify');
            % % % do a multistart fmincon
            optTraj = NaN(predHorizon,numStarts);
            FVAL = NaN(1,numStarts);
            for ii = 1:numStarts
            % % % random trajectory between bounds
            iniGuess = lb + (ub-lb).*rand(size(lb));            
            % % % find optimum trajectory using fmincon
            [optTraj(:,ii),FVAL(ii)] =  fmincon( @(stateTrajectory) ...
                -obj.m_objfForFmincon(sk_k,ck_k,...
                Mk,yk,stateTrajectory,predHorizon),...
                iniGuess,[A3;A4],[b3;b4],[],[],lb(:),ub(:),[],options);
            end
            % % % output
            [val.objFunValFmin,idx] = max(-FVAL);
            val.optStateTrajectoryFmin = optTraj(:,idx)';
            
            fprintf('fmincon optimization complete.\n')
            
            % % % find optimum trajectory using PSO
%             [optTrajPSO,FVALPSO] =  particleswarm( @(stateTrajectory) ...
%                 -obj.m_objfForFmincon(sk_k,ck_k,...
%                 Mk,yk,stateTrajectory,predHorizon),...
%                 predHorizon,lb,ub);
%             % % % output
%             val.optStateTrajectoryPso = optTrajPSO;
%             val.objFunValPso = -FVALPSO;
            
        end
        
        % % % %         objective function for fmincon
        function val = m_objfForFmincon(obj,sk_k,ck_k,Mk,...
                yk,traject,predHorizon)
            
            gpkfPred = obj.m_predictionGPKF(sk_k,ck_k,...
                Mk,yk,traject,predHorizon);
            % % % calculate acquisition function values
            aqFunVal = obj.m_calcAcquisitionFun(gpkfPred.predMeanAtLoc(:),...
                gpkfPred.postVarAtLoc(:),yk);
            aqFunVal = (predHorizon:-1:1).*aqFunVal(:)';
            val = sum(aqFunVal);
            
        end
        
        % % % %         traditional GP regression
        function [predMean,postVar] = traditionalGpRegression(obj,xVisited,yVisited,...
                xPredict,tPredict)
            % % number of points over which we want to acquire predictions
            xPredictNp = size(xPredict,2);
            % % total number of points visited
            xVisitedNp = size(xVisited,2);
            % % form the covariance matrix and mean vector
            [covMat,~] = obj.buildCovMatAndMeanVec(xVisited);
            % % add noise to the covariance
            covMat = covMat + eye(size(covMat))*obj.p_noiseVariance;
            % % preallocate matrices
            predMean = NaN(xPredictNp,1);
            postVar = NaN(xPredictNp,1);
            Kxstar_x = NaN(xPredictNp,xVisitedNp);
            % % prediction points
            xPredictNew = [xPredict;ones(1,xPredictNp)*tPredict];
            % % calculate prediction mean and covariance
            for ii = 1:xPredictNp
                for jj = 1:xVisitedNp
                    Kxstar_x(ii,jj) = obj.m_calcTotCovariance(xPredictNew(:,ii)...
                        ,xVisited(:,jj));
                end
                predMean(ii,1) = Kxstar_x(ii,:)*(covMat\yVisited);
                postVar(ii,1) = obj.m_calcTotCovariance(xPredictNew(:,ii),...
                    xPredictNew(:,ii)) - Kxstar_x(ii,:)*...
                    (covMat\Kxstar_x(ii,:)');
            end
            
        end
        
    end
end

