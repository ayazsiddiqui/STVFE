function [optPts,fMin] = BFGS(rCM,aBooth,bBooth,meanElevation,radius)

    function val = objDist(aBooth,bBooth,meanElevation,radius,s,x,y,z)
        
        val = sqrt((x-radius.*cos((aBooth.*bBooth.^2.*sin(s))./(aBooth.^2.*cos(s).^2+bBooth.^2)).*cos((aBooth.^2.*sin(s.*2.0).*(-1.0./2.0)+bBooth.^2.*meanElevation+aBooth.^2.*meanElevation.*cos(s).^2)./(aBooth.^2.*cos(s).^2+bBooth.^2))).^2+(y+radius.*cos((aBooth.^2.*sin(s.*2.0).*(-1.0./2.0)+bBooth.^2.*meanElevation+aBooth.^2.*meanElevation.*cos(s).^2)./(aBooth.^2.*cos(s).^2+bBooth.^2)).*sin((aBooth.*bBooth.^2.*sin(s))./(aBooth.^2.*cos(s).^2+bBooth.^2))).^2+(z-radius.*sin((aBooth.^2.*sin(s.*2.0).*(-1.0./2.0)+bBooth.^2.*meanElevation+aBooth.^2.*meanElevation.*cos(s).^2)./(aBooth.^2.*cos(s).^2+bBooth.^2))).^2);
    end

x = rCM(1); y = rCM(2); z = rCM(3);
objF = @(s) objDist(aBooth,bBooth,meanElevation,radius,s,x,y,z);

lb = 0;
ub = 2*pi;
maxIter = 100;
bfgsConvergeTol = 1e-3;
bpStep = 0.1;
bpMaxIter = 500;
gradStep = 0.025;
GsConvergeTol = 1e-4;



%% BFGS
interNo = 1;

% Step 1: Starting point
sTest = linspace(0,2*pi,100);
tDis = nan*sTest;
for cc = 1:numel(sTest)
    tDis(cc) = objF(sTest(cc));
end
[~,minIdx] = min(tDis);
x0 = sTest(minIdx);

grad0 = forwardGradient(objF,0,gradStep);
H0 = 1;

% Step 2: Convergence check
while norm(grad0) >= bfgsConvergeTol && interNo < maxIter
    
    % Step 3: Solve set of linear equations
    direction = -H0\grad0(:);
    % Step 4a: Create bounds for alpha_star by using bounding phase
    [alphaLeft,alphaRight] = boundingPhase(objF,x0,direction,bpStep);
    
    % Step 4b: Use golden section to get alpha_star
    alphaStar = goldenSection(objF,x0,direction,alphaLeft,alphaRight,...
        GsConvergeTol);
    
    % Step 5: Get X_new
    x1 = enforceLimits(x0 + alphaStar*direction);
    grad1 = forwardGradient(objF,x1,gradStep);
    
    % Step 6: Update H
    P = x1 - x0;
    Y = grad1 - grad0;
    D = (Y(:)*Y(:)')/(Y(:)'*P(:));
    E = (grad0(:)*grad0(:)')/(grad0(:)'*direction(:));
    
    H0 = H0 + D + E;
    x0 = x1;
    
    if norm(grad1-grad0) < 0.5
        H0 = 1;
    end
    
    grad0 = grad1;
    interNo = interNo+1;
    
end

optPts = x0;
fMin = objF(optPts);

%% secondary functions
%%%% bounding phase
    function [alphaLeft,alphaRight] = boundingPhase(objF,iniPt,direction,stepSize)
        
        maxK = bpMaxIter;
        fVal = NaN(1,maxK);
        alpha = NaN(1,maxK);
        
        if iniPt == lb
            fVal(1) = objF(iniPt);
            fVal(2) = objF(enforceLimits(iniPt + direction*stepSize));
            fVal(3) = objF(enforceLimits(iniPt + 2*direction*stepSize));
        elseif iniPt == ub
            fVal(1) = objF(enforceLimits(iniPt - 2*direction*stepSize));
            fVal(2) = objF(enforceLimits(iniPt - direction*stepSize));
            fVal(3) = objF(iniPt);
        else
            fVal(1) = objF(enforceLimits(iniPt - direction*stepSize));
            fVal(2) = objF(iniPt);
            fVal(3) = objF(enforceLimits(iniPt + direction*stepSize));
        end
        
        % check direction of alha
        if fVal(1) > fVal(2) && fVal(2) > fVal(3)
            cs = 1;
        elseif fVal(1) < fVal(2) && fVal(2) < fVal(3)
            stepSize = -1*stepSize;
            cs = 1;
        else
            cs = 2;
        end
        
        % initialize alpha
        alpha(1:3) = stepSize*[-1;0;1];
        
        switch cs
            case 1
                k = 1;
                while fVal(k+2) <= fVal(k+1) && k<=maxK-2
                    k = k+1;
                    alpha(k+2) = alpha(k+1) + (2^(k-2))*stepSize;
                    fVal(k+2) = objF(enforceLimits(iniPt + direction*alpha(k+2)));
                    
                end
                alphaLeft = alpha(k);
                alphaRight = alpha(k+2);
                
                if k == maxK-2
                    error('bounding phase failed');
                end
            case 2
                alphaLeft = -stepSize;
                alphaRight = stepSize;
        end
    end

%%%% limit checker
    function X = enforceLimits(X)
        belowLim = X<lb;
        X(belowLim) = lb(belowLim);
        aboveLim = X>ub;
        X(aboveLim) = ub(aboveLim);
    end

%%%% golden section
    function alphaStar = goldenSection(objF,iniPt,direction,alphaLeft,alphaRight,...
            convergeTol)
        
        % initial length
        LStart = abs(alphaRight - alphaLeft);
        L(1) = LStart;
        
        % golden ratio
        tau = 0.381966;
        tauI = 1 - tau;
        
        alphaOne = tauI*alphaLeft + tau*alphaRight;
        alphaTwo = tau*alphaLeft + tauI*alphaRight;
        
        f1 = objF(enforceLimits(iniPt + alphaOne*direction));
        f2 = objF(enforceLimits(iniPt + alphaTwo*direction));
        
        k = 1;
        
        while L/LStart > convergeTol && k<1000
            
            if f2 > f1
                
                alphaRight = alphaTwo;
                
                alphaTwo = alphaOne;
                f2 = f1;
                
                alphaOne = tauI*alphaLeft + (tau*alphaRight);
                f1 = objF(enforceLimits(iniPt + (alphaOne*direction)));
                
                L = abs(alphaRight - alphaLeft);
                
                k = k+1;
                
            else
                
                alphaLeft = alphaOne;
                
                alphaOne = alphaTwo;
                f1 = f2;
                
                alphaTwo = (tau*alphaLeft) + tauI*alphaRight;
                f2 = objF(enforceLimits(iniPt + (alphaTwo*direction)));
                
                L = abs(alphaRight - alphaLeft);
                
                k = k+1;
                
            end
        end
        
        alphaStar = (alphaLeft + alphaRight)/2;
        
    end

%%%% forward gradient
    function grad = forwardGradient(objF,iniPt,dX)
        
        grad = NaN*iniPt;
        nEl = numel(iniPt);
        f0 = objF(iniPt);
        fNext = NaN(nEl,1);
        xNext = iniPt;
        
        for ii = 1:nEl
            xNext(ii) = xNext(ii) + dX;
            fNext(ii) = objF(enforceLimits(xNext));
            grad(ii) = (fNext(ii) - f0)/dX;
            xNext = iniPt;
        end
    end

end