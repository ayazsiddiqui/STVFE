function val = windPowerLawMeanFn(z,a,b)
if nargin == 1
    val = z.*(3.17)^0.14;
else
    val = z.*b^a;
end
end

