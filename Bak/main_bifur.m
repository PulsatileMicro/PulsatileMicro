%% Generate input file of arterial bifurcation for 1D program
clear;clc;close all;

VesType='Bifur';

ModelParam=zeros(12,1);
% Measured/Adapted switcher 0: use adapted diameter, 1: use measured diameter
ModelParam(1)=0;
MeasDiam=ModelParam(1);
% Nondimensionalization flag 1: Nondimensionalization
ModelParam(2)=0;
NoDim=ModelParam(2);
% Viscosity update strategy 0: No viscosity update 1: vitro 2: ESL vivo 3: No ESL vivo
ModelParam(3)=0;
% Riemann solver flag 1: Use riemann solver for the BioFlux
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

VesNum=3;
DiamVec=zeros(1,VesNum);
LenVec=zeros(1,VesNum);
VesNumVec=zeros(1,VesNum);
Bifur=zeros(4,VesNum);
BCVal=zeros(1,VesNum);
Bifur=[0 1 1;2 3 2;2 0 0;3 0 0];
BCTypeAll=[];
BCTypeAll=['u' 'B';'B' 'T';'B' 'T'];
% BCTypeAll=['p' 'B';'B' 'T';'B' 'T'];
BCVal(2:3)=0.0;
SegName=zeros(1,VesNum);
BCType=[];
Bifur=Bifur';

DiamRatio_L=0.7583*0.8;
DiamRatio_R=0.7583*0.8;
LenVec=[0.4 0.2 0.2];
DiamVec=[0.02 0.02*DiamRatio_L 0.02*DiamRatio_R];
WallThVec=1e-3;

Vel=100;
VesParam=zeros(20,VesNum);
VesParam(1,:)=LenVec;
VesParam(2,:)=DiamVec;
VesParam(3,:)=WallThVec;         % Wall thickness
VesParam(4,:)=6e5*ones(1,VesNum);  % Young's modulus
VesParam(5,:)=1.26*ones(1,VesNum);  % Alpha
VesParam(6,:)=4e-3;                 % Viscosity
VesParam(7,:)=0.45;                 % Discharge hematocrit
VesParam(8,:)=0;                    % Venous pressure
VesParam(9,:)=3.25E-11;             % Peak of the input flow rate
VesParam(10,:)=0;                   % Segment name
VesParam(11,:)=0;                   % From nodes
VesParam(12,:)=0;                   % To nodes
VesParam(13,:)=Vel;
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
VesParam(21,:)=5;             % Number of history points
VesParam(22,:)=0;             % Taper Rate

uIn=1000;  % mm/s
% Scale the input velocity profile
inFileName='Bifur_IN_1.bcs';
outFileName='Bifur_OUT_1.bcs';
VelProc(dt,inFileName,outFileName,Vel,Vel,Nstep*dt/0.8,1);
% PressureProc(dt,inFileName,outFileName,150.8/10,0,Nstep*dt/0.8,1);
% BoundVel=[Vel*100 1 0]; SinVelProc(dt,inFileName,outFileName,BoundVel,Nstep*dt/0.8,1);
% VelVec=[uIn 1 0]; GaussVelProc(dt,inFileName,outFileName,VelVec,Nstep*dt/0.8,1);

fileName=GenInput(VesType,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);

system('copy Bifur*.* E:\Projects\CPP\1DWin_MKL\1DWin_MKL\');
system('copy Bifur*.* E:\1DWin\Bifur\');
