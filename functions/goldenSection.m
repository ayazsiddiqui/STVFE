function xStar = goldenSection(objF,lb,ub,convergeTol)

xL = min(lb,ub);
xU = max(lb,ub);
% initial length
Ls = xU - xL;
L  = Ls;
% golden ratio
d = (sqrt(5) - 1)/2;

x1 = xL + d*(xU - xL);
x2 = xU - d*(xU - xL);

f1 = objF(x1);
f2 = objF(x2);

k = 1;

while L/Ls > convergeTol && k<26
    
    if f1 < f2
        
        xL = x2;
        
        x2 = x1;
        f2 = f1;
        
        x1 = xL + d*(xU - xL);
        f1 = objF(x1);
        
        L = xU - xL;
        
        k = k+1;
        
    else
        
        xU = x1;
        
        x1 = x2;
        f1 = f2;
        
        x2 = xU - d*(xU - xL);
        f2 = objF(x2);
        
        L = xU - xL;
        
        k = k+1;
        
    end
end

xStar = (xU + xL)/2;

end