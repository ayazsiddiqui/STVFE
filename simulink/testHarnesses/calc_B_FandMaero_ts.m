clear
clc

% load cl and cd values
[AoA_vals,CL_vals,CD_vals] = yClarkCLandCD;
Sref = 10;
rho  = 1.225;
vApp = [-10;0;0];
testAoA = 5;
rollVal = 5*pi/180;

sim('calc_B_FandMaero_th');
