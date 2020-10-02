clear
clc

% simulation parameters
AoA = 5;
init_R = 100;
init_azim = 0;
init_elev = 0;
init_azimRate = 1;
init_elevRate = 0;
gravityEnabled = true;
fluidLoadsEnabled = true;

% kite parameters
kiteMass = 0.7;
refArea  = 0.28;
aeroCent = [0;0;0];
[AoA_vals,CL_vals,CD_vals] = yClarkCLandCD;

% environment
rho = 1.22;
gravAccl = 9.81;
flowSpeed = 10;

%%
tiledlayout(2,1)

nexttile
plot(AoA_vals,CL_vals);
grid on; hold on;
xlabel('AoA [deg]'); ylabel('$C_{L}$');

nexttile
plot(AoA_vals,CD_vals);
grid on; hold on;
xlabel('AoA [deg]'); ylabel('$C_{D}$');
