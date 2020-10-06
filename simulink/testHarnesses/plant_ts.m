clear
clc

% path parameters
pathWidth = 50;
pathHeight = 10;
meanElevation = 30;

pathBasis = calcPathBasisVector(pathWidth,pathHeight,meanElevation);
pathBasis = timeseries(repmat(pathBasis(:),1,2),[0 1e6]);

% simulation parameters
AoA = 5;
init_R = 100;
init_azim = 0;
init_elev = meanElevation*pi/180;


init_Speed = 10;
init_Heading = 210;
init_azimRate = init_Speed*cosd(init_Heading)/init_R;
init_elevRate = init_Speed*sind(init_Heading)/init_R;

% kite parameters
kiteMass = 5;
refArea  = 0.28;
aeroCent = [0;0;0];
[AoA_vals,CL_vals,CD_vals] = yClarkCLandCD;
% environment
rho = 1.22;
gravAccl = 9.81;
flowSpeed = 10;
rollAng = 30*pi/180;


%% simulation
simTime = 10;
out = sim('plant_th');

%% animate
tsc = parseLogsout(out,'resample',1);
pp = PLOT.plottingSuite;

figure
pp.makeAnimation(tsc.G_rCM,tsc.G_bdy,tsc.pathBasis,'G_vCM',tsc.G_vCM);


%%
% tiledlayout(2,1)
% 
% nexttile
% plot(AoA_vals,CL_vals);
% grid on; hold on;
% xlabel('AoA [deg]'); ylabel('$C_{L}$');
% 
% nexttile
% plot(AoA_vals,CD_vals);
% grid on; hold on;
% xlabel('AoA [deg]'); ylabel('$C_{D}$');
