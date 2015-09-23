clear;clc;close all;

ModelParam=zeros(12,1);
% Measured/Adapted switcher 0: use adapted diameter, 1: use measured diameter
ModelParam(1)=0;
MeasDiam=ModelParam(1);
% Nondimensionalization flag 1: Nondimensionalization
ModelParam(2)=0;
NoDim=ModelParam(2);
% Viscosity update strategy 0: No viscosity update 1: No ESL considered 2: ESL considered
ModelParam(3)=0;
% Riemann solver flag 1: Use Riemann solver for the BioFlux
ModelParam(4)=1;
% Binary output flag 0: Output text format 1: Output binary format
ModelParam(5)=0;
% Time step
ModelParam(6)=1e-3;
dt=ModelParam(6);
% Number of steps
ModelParam(7)=6*0.8/dt;
Nstep=ModelParam(7);
% Viscoelastic wall 0: No viscous wall 1: New viscous model 2: Old viscous model
ModelParam(9)=0;
% Show CFL flag 0: Hide CFL number 1: Show CFL flag
ModelParam(10)=0;
% Show Step Lapse flag 0: Hide Step Lapse 1: Show Step Lapse
ModelParam(11)=0;
% ODE solver flag 0: Explicit 1:Implicit
ModelParam(12)=1;

% order指从起点到毛细血管层的阶数，总阶数为2*order-1
% order=2时即为bifur
order=5;
curOrder=0;
cnt=0;
lastNum=0;
VesNum=2^order-1;

% Parameters of start vessel
% startDiam=4e-5/0.71;
startDiam=4e-5;
ratio_L=60;

% ratio_L=20;
% ratio_L=50*0.85/0.77*0.74/0.70;
% ratio_D_L=0.7583;
% ratio_D_R=0.7583;
ratio_D_L=0.89;
ratio_D_R=0.89;

DiamVec=zeros(1,VesNum);
LenVec=zeros(1,VesNum);
VesNumVec=zeros(1,VesNum);
Bifur=zeros(4,VesNum);
BCTypeAll=[];
SegName=zeros(1,VesNum);
BCVal=zeros(1,VesNum);

for i=1:VesNum
  BCType=[];
  VesNumVec(i)=i;
  SegName(i)=i;
  if mod(i,2^(curOrder))==0
    curOrder=curOrder+1;
  end
  if curOrder<=order
    if i==1
      DiamVec(1)=startDiam;
      LenVec(1)=startDiam*ratio_L;
      Bifur(1:2,1)=[0 2];
    elseif mod(i,2)==0
      DiamVec(i)=DiamVec(i/2)*ratio_D_L;
      LenVec(i)=DiamVec(i)*ratio_L;
      Bifur(1,i)=i/2;
      Bifur(2,i)=i+1;
      Bifur(3,i/2)=i;
      Bifur(4,i/2)=i+1;
    else
      DiamVec(i)=DiamVec((i-1)/2)*ratio_D_R;
      LenVec(i)=DiamVec(i)*ratio_L;
      Bifur(1,i)=(i-1)/2;
      Bifur(2,i)=i-1;
    end
    if curOrder==order
      BCType=[BCType 'B' 'T'];
      BCVal(i)=0.8;
    elseif i == 1
      BCType=[BCType 'u' 'B'];
%       BCType=[BCType 'p' 'B'];
    else
      BCType=[BCType 'B' 'B'];
    end
  end
  BCTypeAll=[BCTypeAll;BCType];
end
Bifur=Bifur';
LenVec=[0.00240000000000000,0.00213600000000000,0.00213600000000000,0.00190104000000000,0.00190104000000000,0.00190104000000000,0.00190104000000000,0.00169192560000000,0.00169192560000000,0.00169192560000000,0.00169192560000000,0.00169192560000000,0.00169192560000000,0.00169192560000000,0.00169192560000000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000,0.00150581378400000;];
% LenVec=[0.4 0.2 0.2]*1e-2;

Vel=100;
VesParam=zeros(22,VesNum);
VesParam(1,:)=LenVec*5;
VesParam(2,:)=DiamVec;
VesParam(3,:)=startDiam*0.1;         % Wall thickness
VesParam(4,:)=3.5e5*ones(1,VesNum);  % Young's modulus
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
VesParam(16,:)=3;                   % L order
VesParam(17,:)=3;                   % q order
if NoDim
    VesParam(18,:)=1e-3;                % scale_lamda
    VesParam(19,:)=1e-1;                % scale_u0
    VesParam(20,:)=1e-4;                % scale_r0
else
    VesParam(18,:)=1;                % scale_lamda
    VesParam(19,:)=1;                % scale_u0
    VesParam(20,:)=1;                % scale_r0
end
VesParam(21,:)=3;                    % Number of history points
VesParam(22,:)=0;                    % Taper Rate

VesType='Tree';
inFileName='Tree_IN_1.bcs';
outFileName='Tree_OUT_1.bcs';
% Scale the input velocity profile
VelProc(dt,inFileName,outFileName,Vel,Vel,Nstep*dt/0.8,1);
% VelVec=[1e4 1 0]; SinVelProc(dt,inFileName,outFileName,VelVec,Nstep*dt/0.8,1);
% PressureProc(dt,inFileName,outFileName,Vel,Vel,Nstep*dt/0.8,1);

% fileName=GenSymNet(VesType, VesParam, BCTypeAll, Bifur', dt, Nstep);
fileName=GenInput(VesType,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);

system('copy Tree*.* E:\Projects\CPP\1DWin_MKL\1DWin_MKL\');
system('copy Tree*.* E:\1DWin\Tree\');
