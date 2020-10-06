function val = calcPathBasisVector(pathWidth,pathHeight,pathMeanElevation)

% get bBooth
% local variables
w = pathWidth*pi/180;
h = pathHeight*pi/180;
% output
val(1) = 0.5*w;
val(2) = (1/(2*sqrt(2)))*sqrt(-w^2+sqrt((h^2*(4+h^2)*w^4))/(h^2));
val(3) = pathMeanElevation*pi/180;

end
