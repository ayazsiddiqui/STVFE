function val = getPathCoords(aBooth,bBooth,meanElevation,radius,s)

s = 2*pi - s;
phi = aBooth*sin(s)/(1 + (aBooth*cos(s)/bBooth)^2);
beta = ((aBooth/bBooth)^2*sin(s)*cos(s)/(1 + (aBooth*cos(s)/bBooth)^2)) + meanElevation;

val = radius*[cos(phi)*cos(beta);sin(phi)*cos(beta);sin(beta)];

end