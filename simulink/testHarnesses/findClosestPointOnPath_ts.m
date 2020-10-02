clear
clc

% parameters
G_rCM = [75;45;100];
aBooth = 0.4363;
bBooth = 0.9999;
meanElev = 30*pi/180;
radius = 150;
s0 = 0.25*2*pi;
L1 = 200;

options = optimoptions('fmincon','algorithm','sqp');
sBest = fmincon( @(s) ...
    distBetweenPts(G_rCM,s,aBooth,bBooth,meanElev,radius),s0,[],[],[],[],...
    s0,s0 + 0.25*2*pi,[],options);

rClosest = getPathCoords(aBooth,bBooth,meanElev,radius,sBest);

[sIntersect,minVal] = fmincon(@(ds)...
    calcTargetDistance(aBooth,bBooth,meanElev,radius,s0,L1,ds),sBest,...
    [],[],[],[],0,[],[],options);

rIntersect = getPathCoords(aBooth,bBooth,meanElev,radius,sIntersect);

sTest = linspace(0,2*pi,201);

pathXYZ = getPathCoords(aBooth,bBooth,meanElev,radius,sTest);
plot3(pathXYZ(1,:),pathXYZ(2,:),pathXYZ(3,:),'k-');
grid on;hold on;
view(120,30);
plot3(rClosest(1),rClosest(2),rClosest(3),'r*');
plot3(G_rCM(1),G_rCM(2),G_rCM(3),'b*');
plot3(rIntersect(1),rIntersect(2),rIntersect(3),'m*');
axis equal;