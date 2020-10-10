clear
clc

objF = @(x) -4.*sin(x).*(1 + cos(x));
% objF = @(x) (x.^2 + 2);
lb = -1;
ub = 3;
convergeTol = 1e-5;

xS = goldenSection(objF,lb,ub,convergeTol);


xT = linspace(lb,ub,101);
plot(xT,objF(xT));
hold on
plot(xS,objF(xS),'r*');

