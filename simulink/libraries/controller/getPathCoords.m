function val = getPathCoords(aBooth,bBooth,meanElevation,radius,s)

s = 2*pi - s;
phi = (aBooth.*sin(s))./(aBooth.^2.*1.0./bBooth.^2.*cos(s).^2+1.0);
beta = meanElevation+(aBooth.^2.*1.0./bBooth.^2.*cos(s).*sin(s))./(aBooth.^2.*1.0./bBooth.^2.*cos(s).^2+1.0);

val = radius.*[cos(phi).*cos(beta);sin(phi).*cos(beta);sin(beta)];

end