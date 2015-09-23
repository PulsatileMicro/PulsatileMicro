%% Main function
clear;clc;close all;

VesType='Net_122';    % Network name
DatFile='test.DAT';   % Network data file name
PrnFile='test.prn';

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
ModelParam(7)=4*Period/dt;
Nstep=ModelParam(7);
% Viscoelastic wall 0: No viscous wall 1: New viscous model 2: Old viscous model
ModelParam(9)=0;
% Show CFL flag 0: Hide CFL number 1: Show CFL flag
ModelParam(10)=0;
% Show Step Lapse flag 0: Hide Step Lapse 1: Show Step Lapse
ModelParam(11)=1;
% ODE solver flag 0: Explicit 1:Implicit
ModelParam(12)=1;

% Read data from input files
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
SegName=DataArray(:,1);
VesCategory=DataArray(:,2);
From=DataArray(:,3);
To=DataArray(:,4);
Len=DataArray(:,6);
if MeasDiam==0
  Diam=FuncPara(:,3); % Read adapted diameter
else
  % 122 Network doesn't have measured diameter
end
Visc=FuncPara(:,9);
WallTh=FuncPara(:,20);
Hd=FuncPara(:,10);
if MeasDiam==0
  Vel=FuncPara(:,7);
else
  % 122 Network doesn't have measured diameter
end

% Get number of vessels
VesNum=length(SegName);
BCTypeAll=[];   % Boundary condition type vector
Bifur=zeros(VesNum,4);
BCVal=zeros(VesNum,1);
% Set the boundary type, bifurcation condition and value
for i=1:VesNum
  BCType=[];
  % Determine the inlet type
  inInd1=find(From(i)==From);
  inInd2=find(From(i)==To);
  if length(inInd1)==2
    BCType=[BCType 'B'];
    inInd1(inInd1==i)=[];
    Bifur(i,1:2)=[inInd2 inInd1];
  elseif length(inInd2)==2
    BCType=[BCType 'C'];
    Bifur(i,1:2)=inInd2;
  elseif length(inInd2)==0
    %         BCType=[BCType 'u'];
    BCType=[BCType 'q'];
    Bifur(i,1:2)=[0 2];
  else
    BCType=[BCType 'J'];
    Bifur(i,1:2)=inInd2;
  end
  % Determine the outlet type
  outInd1=find(To(i)==From);
  outInd2=find(To(i)==To);
  if length(outInd1)==2
    BCType=[BCType 'B'];
    Bifur(i,3:4)=outInd1;
  elseif length(outInd2)==2
    BCType=[BCType 'C'];
    outInd2(outInd2==i)=[];
    Bifur(i,3:4)=[outInd2 outInd1];
  elseif length(outInd1)==0
    % Only ONE output node in this network
    Bifur(i,3:4)=[0 0];
    %         BCType=[BCType 'R'];
    %         BCVal(i)=2.7531e+14;
    BCType=[BCType 'T'];
    BCVal(i)=0;
  else
    BCType=[BCType 'J'];
    Bifur(i,3:4)=outInd1;
  end
  BCTypeAll=[BCTypeAll;BCType];
end

% Map the vessel order to the real vessel number in the original dataset
% It facilitates the validation by comparing with the figure of the network
BifurSeg=zeros(VesNum,4);
for i=1:VesNum
  for j=1:4
    if Bifur(i,j)~=0
      BifurSeg(i,j)=SegName(Bifur(i,j));
    else
      BifurSeg(i,j)=0;
    end
  end
end

VesParam=zeros(19,VesNum);

% kg,m,s
VesParam(1,:)=Len'*1e-6;            % Length
VesParam(2,:)=Diam'*1e-6;           % Diameter
VesParam(3,:)=WallTh'*1e-6;         % Wall thickness
% VesParam(3,:)=mean(WallTh)*1e-6;         % Wall thickness
for i=1:VesNum
  if VesCategory(i)==1
    VesParam(4,i)=3.5e5;
  elseif VesCategory(i)==3
    VesParam(4,i)=3.88e5;
  else
    VesParam(4,i)=3.7e5;
  end
end
VesParam(5,:)=1.26*ones(1,VesNum);  % Alpha
VesParam(6,:)=Visc'*1e-3;           % Viscosity
% VesParam(6,:)=Visc'*1e-4;           % Viscosity
VesParam(7,:)=Hd;                   % Discharge hematocrit
VesParam(8,:)=0;                    % Venous pressure
VesParam(9,:)=0;                    % Peak of the input flow rate
VesParam(10,:)=SegName;             % Segment name
VesParam(11,:)=From;                % From nodes
VesParam(12,:)=To;                  % To nodes
VesParam(13,:)=Vel;
if ModelParam(9)==1
  VesParam(14,:)=Eval_Gamma(VesParam(2,:),1e4*ones(1,VesNum),VesParam(3,:));  % Visc part of viscoelasticity
elseif ModelParam(9)==2
  VesParam(15,:)=Eval_GammaII(VesParam(2,:),1e4*ones(1,VesNum));  % Visc part of viscoelasticity
else
  VesParam(14,:)=0;
  VesParam(15,:)=0;
end
VesParam(16,:)=3;                   % q
VesParam(17,:)=3;                   % L
if NoDim==1
  VesParam(18,:)=1e-3;                % scale_lamda
  VesParam(19,:)=1e-1;                % scale_u0
  VesParam(20,:)=1e-4;                % scale_r0
  %     VesParam(18,:)=1050;                % scale_lamda
  %     VesParam(19,:)=1;                % scale_u0
  %     VesParam(20,:)=1;                % scale_r0
else
  VesParam(18,:)=1;
  VesParam(19,:)=1;
  VesParam(20,:)=1;
end

inFileName='Net_122_IN_1.bcs';
outFileName='Net_122_OUT_84.bcs';
% Scale the input velocity profile
% MeanBottomRatioA=VelProc(dt,inFileName,outFileName,Vel(1),Vel(84),Nstep*dt/0.8,VesParam(18,1));
MeanBottomRatioA=FlowProc(dt,inFileName,outFileName,Vel(1),Vel(84),Nstep*dt/0.8,Diam,VesParam(18,1),VesParam(19,1));
ModelParam(8)=MeanBottomRatioA;

% Generate .in file of the simulation
fileName=GenInput(VesType,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);
system('copy Net_122*.* E:\Projects\CPP\1DWin_MKL\1DWin_MKL\');
system('copy Net_122*.* E:\1DWin\122\');
