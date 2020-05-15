classdef GPKF
    %GPKF Summary of this class goes here
    %   Detailed explanation goes here
    
    %% properties
    properties
        noSpatialIps
    end
    
    %% methods
    methods
        
        % % % %         constructor
        function obj = GPKF(noSpatialIps)
            obj.noSpatialIps = noSpatialIps;
        end
        
        % % % %         spatial kernel: squared exponential
        function val = meanFunction(obj,x)
            % % zero mean function
            val = 0*(x'*x)*obj.noSpatialIps;
        end
        
        % % % %         spatial kernel: squared exponential
        function val = spatialKernel(obj,s1,s2,covAmp,lengthScales)
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
                    covMat(ii,jj) = obj.spatialKernel(x(:,ii),x(:,jj),...
                        covAmp,lengthScales);
                end
            end
            % form the total covariance matrix
            covMat = covMat + triu(covMat,1)';
        end
        
        % % % %         temporal kernel: exponential
        function val = temporalKernel(obj,t1,t2,timeScale)
            % % covariance equation
            val = 1*exp(-0.5*((t1-t2)^2)/timeScale^2);
        end
        
        % % % %         calculate covaraince as a product of the two covariances
        function val = calcTotCovariance(obj,x1,x2,hyperParams)
            % % covariance amplitude or variance of latent function
            covAmp = hyperParams(1);
            % % length scales for spatial covariance
            lenScale = hyperParams(2:obj.noSpatialIps+1);
            % % time scale
            timeScale = hyperParams(obj.noSpatialIps+2);
            % % k = k_s*k_t
            val = obj.spatialKernel(x1(1:end-1),x2(1:end-1),covAmp,lenScale)...
                *obj.temporalKernel(x1(end),x2(end),timeScale);
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
        
        % % % %         exponential kernel GPKF initialization
        function val = exponentialGpkfInitialize(obj,xDomain,timeScale,timeStep)
            % % total number of points in the entire domain of interest
            xDomainNP = size(xDomain,2);
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
        
        % % % %         squared exponential kernel GPKF initialization
        function val = squaredExponentialGpkfInitialize(obj,xDomain,timeScale,...
                timeStep,N)
            % % total number of points in the entire domain of interest
            xDomainNP = size(xDomain,2);
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
            coEffs = obj.removeEPS(coEffs,6);
            % % normalize them by dividing by the highest degree
            % % coefficient
            coEffs = coEffs./coEffs(end);
            % % form the F, G, and H matrices as per Carron Eqn. (8)
            F = [zeros(N-1,1) eye(N-1); -coEffs(1:end-1)];
            G = [zeros(N-1,1);1];
            % % calculate the numerator
            b0 = sqrt((timeScale^2)*factorial(N)*((2/(timeScale^2))^N)...
                *sqrt(pi*2*timeScale^2));
            H = [b0 zeros(1,N-1)];
            sigma0 = obj.removeEPS(lyap(F,G*G'),6);
            % % calculate the discretized values
            syms tau
            % % use cayley hamilton theorem to calcualte e^Ft
            eFt = obj.cayleyHamilton(F);
            % % calculate Fbar using the above expression
            Fbar = obj.removeEPS(subs(eFt,tau,timeStep),6);
            % % evaluate Qbar, very computationally expensive
            Qsym = eFt*(G*G')*eFt';
            Qint = NaN(N);
            for ii = 1:N^2
                fun = matlabFunction(Qsym(ii));
                Qint(ii) = integral(fun,0,timeStep);
            end
            % % remove numbers lower than eps
            Qbar = obj.removeEPS(Qint,6);
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
        
        % % % %         Kalman estimation as per jp Algorithm 1
        function [F_t,sigF_t,skp1_kp1,ckp1_kp1] = ...
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
            skp1_kp1 = skp1_k + Lkp1*(yk - Cmat*skp1_k); % Eqn (6c)
            ckp1_kp1 = ckp1_k - Lkp1*Cmat*ckp1_k; % Eqn (6d)
            % % process estimate and covariance as per Todescato algortihm 1
            F_t = Ks_12*Hmat*skp1_kp1; % Eqn. (13)
            sigF_t = Ks_12*Hmat*ckp1_kp1*Hmat'*Ks_12;
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
            h0 = obj.temporalKernel(0,0,timeScale);
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
                sigmaX(ii,jj) = h0*obj.spatialKernel(xPredict(:,ii),xDomain(:,jj)...
                    ,covAmp,lengthScales);
                end
                Vx(ii,1) = h0*obj.spatialKernel(xPredict(:,ii),xPredict(:,ii)...
                    ,covAmp,lengthScales);
                % % predicted mean as per Todescato Eqn. (17)
                predMean(ii,1) = sigmaX(ii,:)*(Vf\F_t);
                % % posterior variance as per Todescato Eqn. (18)
                postVar(ii,:) = Vx(ii,1) - sigmaX(ii,:)*(Vf\(Vf - sigF_t))*...
                    (Vf\sigmaX(ii,:)');
            end
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

