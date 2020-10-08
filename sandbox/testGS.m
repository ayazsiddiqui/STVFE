clear
clc

objF = @(x) -sin(x);
rCM = [100;0;150];

aBooth = 0.5;
bBooth = 1;
meanElevation = 30*pi/180;
radius = 150;
s0 = 0;
convergeTol = 1e-5;

