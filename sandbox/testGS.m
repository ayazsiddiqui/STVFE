clear
clc

objF = @(x) -sin(x);

iniPt = 2*pi;
alphaLeft = iniPt;
alphaRight = iniPt + 2*pi;
convergeTol = 1e-6;
direction = -cos(iniPt);

alphaStar = goldenSection(objF,iniPt,direction,alphaLeft,alphaRight,...
            convergeTol);


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
        
        f1 = objF(iniPt + alphaOne*direction);
        f2 = objF(iniPt + alphaTwo*direction);
        
        k = 1;
        
        while L/LStart > convergeTol && k<1000
            
            if f2 > f1
                
                alphaRight = alphaTwo;
                
                alphaTwo = alphaOne;
                f2 = f1;
                
                alphaOne = tauI*alphaLeft + (tau*alphaRight);
                f1 = objF(iniPt + (alphaOne*direction));
                
                L = abs(alphaRight - alphaLeft);
                
                k = k+1;
                
            else
                
                alphaLeft = alphaOne;
                
                alphaOne = alphaTwo;
                f1 = f2;
                
                alphaTwo = (tau*alphaLeft) + tauI*alphaRight;
                f2 = objF(iniPt + (alphaTwo*direction));
                
                L = abs(alphaRight - alphaLeft);
                
                k = k+1;
                
            end
        end
        
        alphaStar = (alphaLeft + alphaRight)/2;
        
    end