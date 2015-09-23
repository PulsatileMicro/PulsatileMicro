%% Generate input file of single artery for 1D program
clear;clc;close all;

% Network name
VesType='Artery';
VesNum=1;
SegName=1;

ModelParam=zeros(12,1);
% Measured/Adapted switcher 0: use adapted diameter, 1: use measured diameter
ModelParam(1)=0;
MeasDiam=ModelParam(1);
% Nondimensionalization flag 1: Nondimensionalization
ModelParam(2)=0;
NoDim=ModelParam(2);
% Viscosity update strategy 0: No viscosity update 1: vitro 2: ESL vivo 3: No ESL vivo
ModelParam(3)=0;
% Riemann solver flag 1: Use Riemann solver for the BioFlux
ModelParam(4)=1;
% Binary output flag 0: Output text format 1: Output binary format
ModelParam(5)=0;
% Time step
ModelParam(6)=1e-3;
dt=ModelParam(6);
% Number of steps
Period=0.8;
% Period=0.325;
ModelParam(7)=6*Period/dt;
Nstep=ModelParam(7);
% Viscoelastic wall 0: No viscous wall 1: New viscous model 2: Old viscous model
ModelParam(9)=0;
% Show CFL flag 0: Hide CFL number 1: Show CFL flag
ModelParam(10)=0;
% Show Step Lapse flag 0: Hide Step Lapse 1: Show Step Lapse
ModelParam(11)=0;
% ODE solver flag 0: Explicit 1:Implicit
ModelParam(12)=1;

VesParam=zeros(21,VesNum);
% data from 
% Wang, J. J. and K. H. Parker. Wave propagation in a model of the arterial circulation. J Biomech. 37(4):457-470, 2004.
% R. femoral
% VesParam(1,:)=0.443;          % Length
% VesParam(2,:)=0.00762;        % Diameter
% VesParam(3,:)=0.0005;         % Wall thickness
% VesParam(4,:)=0.8e6;          % Young's modulus
% Test data
VesParam(1,:)=0.4;          % Length
VesParam(2,:)=0.02;         % Diameter
VesParam(3,:)=0.001;        % Wall thickness
VesParam(4,:)=0.6e6;          % Young's modulus
VesParam(5,:)=1.1;            % Alpha
VesParam(6,:)=4e-3;           % Viscosity

% VesParam(1,:)=0.02;
% VesParam(2,:)=4e-5;         % Diameter
% VesParam(3,:)=4e-6;        % Wall thickness
% VesParam(4,:)=3.5e5;          % Young's modulus
% VesParam(5,:)=1.26;            % Alpha
% VesParam(6,:)=2e-3;           % Viscosity

% VesParam(1,:)=0.4/20;
% VesParam(2,:)=0.02/500;
% VesParam(3,:)=0.001/500;
% VesParam(4,:)=3.5e5;          % Young's modulus
% VesParam(5,:)=1.26;            % Alpha
% VesParam(6,:)=2e-3;           % Viscosity

VesParam(7,:)=0.45;           % Discharge hematocrit
VesParam(8,:)=0;              % Venous pressure
VesParam(9,:)=312.7E-6;       % Peak of the input flow rate
VesParam(13,:)=20;
if ModelParam(9)==1
  VesParam(14,:)=Eval_Gamma(VesParam(2,:),1e4*ones(1,VesNum),VesParam(3,:));  % Visc part of viscoelasticity
elseif ModelParam(9)==2
  VesParam(15,:)=Eval_GammaII(VesParam(2,:),1e4*ones(1,VesNum));  % Visc part of viscoelasticity
else
  VesParam(14,:)=0;
  VesParam(15,:)=0;
end
VesParam(16,:)=20;            % q
VesParam(17,:)=20;            % L
if NoDim
    VesParam(18,:)=1e-3;                % scale_lamda
    VesParam(19,:)=1e-1;                % scale_u0
    VesParam(20,:)=1e-4;                % scale_r0
else
    VesParam(18,:)=1;                % scale_lamda
    VesParam(19,:)=1;                % scale_u0
    VesParam(20,:)=1;                % scale_r0
end
VesParam(21,:)=5;         % Number of history points
VesParam(22,:)=0;      % Taper Rate

uIn=0.5;  % mm/s
pIn=10;   % mmHg
qIn=100;    % ml/s
BCTypeAll=[];
BCType=[];
BCType=['u' 'T'];
Bifur=[0 2];
BCVal=0.5;
% BCType=['p' 'T'];
% Bifur=[0 2];
% BCVal=0.2;
BCTypeAll=[BCTypeAll;BCType];

inFileName='Artery_IN.bcs';
outFileName='Artery_OUT.bcs';
% Scale the input velocity profile
% PressureProc(dt,inFileName,outFileName,150.8,0,Nstep*dt/0.8,1);   % Input Unit mmHg
VelProc(dt,inFileName,outFileName,uIn,0,Nstep*dt/0.8,1);
% FlowProc(dt,inFileName,outFileName,qIn,0,Nstep*dt/0.8,VesParam(2,:)*1e6,1,1);
% VelVec=[2e8/1e4 1 0]; SinVelProc(dt,inFileName,outFileName,VelVec,Nstep*dt/0.8,1);
% VelVec=[uIn 1 0]; GaussVelProc(dt,inFileName,outFileName,VelVec,Nstep*dt/0.8,1);

% Generate .in file of the simulation
fileName=GenInput(VesType,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);
system('copy Artery*.* E:\Projects\CPP\1DWin_MKL\1DWin_MKL\');
system('copy Artery*.* E:\1DWin\Artery\');
