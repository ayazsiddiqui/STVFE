% clear
% clc
% format compact

cd(fileparts(mfilename('fullpath')));


VEHICLE               = "vehicleLE";
PLANT                 = "plantDOE";
SIXDOFDYNAMICS        = "sixDoFDynamicsCoupled";

%% Essential Values
vhcl = OCT.vehicle;

vhcl.setFluidDensity(1000,'kg/m^3')
vhcl.setNumTethers(1,'');
vhcl.setBuoyFactor(1.0,'');
vhcl.setFluidCoeffsFileName('ayazAirborneVhclCoeffs_v1','');

%% Volumes and Inertia
vhcl.setVolume(2.3385,'m^3') %From CAD
Ixx = 9.143165146555E+09*1e-6;
Iyy = 5.942254391763E+09*1e-6;
Izz = 1.5005672119429E+10*1e-6;
Ixy = 5.742*1e-6;
Ixz = -120250375.384*1e-6;
Iyz = 0;

vhcl.setInertia_CM([Ixx -Ixy -Ixz;...
                    -Ixy Iyy -Iyz;...
                    -Ixz -Iyz Izz],'kg*m^2');
                
%% Important Points
vhcl.setRB_LE([0 0 0],'m');
vhcl.setRCM_LE([0.831 0 0.0224],'m');
vhcl.setRBridle_LE(vhcl.rCM_LE.Value + [0;0;0],'m');
vhcl.setRCentOfBuoy_LE([0.831 0 0.0224],'m');% [.492,0,.0682] from CAD

%% Added Mass/Damping (defaults to zeros)
vhcl.setMa6x6_LE(-0*[134         0           0           0          14           0;...
                   0      1591           0        -986           0        5166;...
                   0         0        9265           0       -8202           0;...
                   0      -986           0       67412           0       -4608;...
                  14         0       -8202           0       20209           0;...
                   0      5166           0       -4608           0       23320;],'');
% vhcl.setD6x6_B([],'');

%% Control Surfaces
vhcl.setAllMaxCtrlDef(30,'deg');
vhcl.setAllMinCtrlDef(-30,'deg');
vhcl.setAllMaxCtrlDefSpeed(60,'deg/s');

% super high clMax so that none of the surfaces stall
clMax = 5;

%% Wing
vhcl.setWingRootChord(1.392,'m');
vhcl.setWingAR(10,''); %Span 10, hspan 5
vhcl.setWingTR(0.315,'');
vhcl.setWingSweep(0,'deg');
vhcl.setWingDihedral(0,'deg');
vhcl.setWingIncidence(0,'deg');
vhcl.setWingAirfoil('n63412.dat','');
vhcl.setWingClMin(-clMax,'');
vhcl.setWingClMax(clMax,'');

%% H-stab and V-stab
vhcl.hStab.setRSurfLE_WingLEBdy([4.071;0;0],'m');
vhcl.hStab.setNumTraps(2,'');
vhcl.hStab.setRootChord(.6429,'m');
vhcl.hStab.setTR(.817,'');
vhcl.hStab.setHalfSpanGivenAR(4.93,''); %Span 4, hspan 2
vhcl.hStab.setSweep(4.24,'deg');
vhcl.hStab.setIncidence(0,'deg');
vhcl.hStab.setAirfoil('NACA0012','');
vhcl.hStab.setClMin(-clMax,'');
vhcl.hStab.setClMax(clMax,'');

vhcl.vStab.setRSurfLE_WingLEBdy([3.513;0;0],'m');
vhcl.vStab.setRootChord(1.2,'m');
vhcl.vStab.setHalfSpan(1.2,'m');
vhcl.vStab.setTR(.446,'');
vhcl.vStab.setSweep(40,'deg');
vhcl.vStab.setAirfoil('NACA0015','');
vhcl.vStab.setClMin(-clMax,'');
vhcl.vStab.setClMax(clMax,'');

%% Fuselage (could use more realistic numbers)
vhcl.fuse.setDiameter(0.2,'m');
vhcl.fuse.setEndDragCoeff(0.0,'');
vhcl.fuse.setSideDragCoeff(0.0,'');
vhcl.fuse.setRNose_LE([-2.785;0;0],'m');
vhclEndPoint = max(vhcl.hStab.rSurfLE_WingLEBdy.Value(1)+vhcl.hStab.rootChord.Value,...
    vhcl.vStab.rSurfLE_WingLEBdy.Value(1)+vhcl.vStab.rootChord.Value);
vhcl.fuse.setREnd_LE([vhclEndPoint;0;0],'m');

%% Turbines
vhcl.setNumTurbines(2,'');
vhcl.build('TurbClass','turb');
% port rotor
vhcl.turb1.setMass(0,'kg')
vhcl.turb1.setDiameter(0,'m')
vhcl.turb1.setAxisUnitVec([1;0;0],'')
vhcl.turb1.setAttachPtVec(vhcl.portWing.outlinePtsBdy.Value(:,2),'m')
vhcl.turb1.setPowerCoeff(.5,'')
vhcl.turb1.setAxalInductionFactor(1.5,'')
vhcl.turb1.setTipSpeedRatio(6,'')
% starboard rotor
vhcl.turb2.setMass(0,'kg')
vhcl.turb2.setDiameter(0,'m')
vhcl.turb2.setAxisUnitVec([-1;0;0],'')
vhcl.turb2.setAttachPtVec(vhcl.stbdWing.outlinePtsBdy.Value(:,2),'m')
vhcl.turb2.setPowerCoeff(.5,'')
vhcl.turb2.setAxalInductionFactor(1.5,'')
vhcl.turb2.setTipSpeedRatio(6,'')
    
%% load/generate fluid dynamic datan
vhcl.calcFluidDynamicCoefffs

%% scale
densityScaleFactor = 35/vhcl.mass.Value;
vhcl.scale(1,densityScaleFactor);

%% save file in its respective directory
saveBuildFile('vhcl',mfilename,'variant',["VEHICLE","PLANT","SIXDOFDYNAMICS"]);



