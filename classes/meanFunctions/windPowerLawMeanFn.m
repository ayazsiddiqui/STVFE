function val = windPowerLawMeanFn(z,a,b)
if nargin == 1
    val = 3.77*z.^0.14;
else
    val = a*z.^b;
end
end

