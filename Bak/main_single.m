%% Generate input file for 1D program
clear;clc;close all;

% Network name
VesType='Single';
% Network data file name
DatFile='Single_2.DAT';
PrnFile='Single_2.prn';
ModelParam=zeros(12,1);
% Measured/Adapted switcher 0: use adapted diameter, 1: use measured diameter
ModelParam(1)=0;
MeasDiam=ModelParam(1);
% Nondimensionalization flag 1: Nondimensionalization
ModelParam(2)=0;
NoDim=ModelParam(2);
% Viscosity update strategy 0: No viscosity update 1: No ESL considered 2: ESL considered
ModelParam(3)=0;
% Riemann solver flag 1: Use riemann solver for the BioFlux
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

% Read data from input files
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
SegName=DataArray(:,1);
VesselType=DataArray(:,2);
From=DataArray(:,3);
To=DataArray(:,4);
Len=DataArray(:,6);
Diam=DataArray(:,5);    % The measured diameter
Visc=FuncPara(:,9);
WallTh=FuncPara(:,20);
% Hd=FuncPara(:,10);
Hd=DataArray(:,7);
% Vel=FuncPara(:,7);
Vel=DataArray(:,8);
MeanP=FuncPara(:,4);
DeltaP=FuncPara(:,24);
Flow=FuncPara(:,8);
VesCatogory=FuncPara(:,19);

% Get number of vessels
VesNum=length(SegName);

BCTypeAll=[];
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
    BCType=[BCType 'u'];
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
    BCType=[BCType 'T'];
    %         BCType=[BCType 'R'];
    Bifur(i,3:4)=[0 0];
    BCVal(i)=0.0;
    %         BCVal(i)=(MeanP(i)-DeltaP(i)/2)/Flow(i)*133*60*1e12;
  else
    BCType=[BCType 'J'];
    Bifur(i,3:4)=outInd1;
  end
  BCTypeAll=[BCTypeAll;BCType];
end

VesParam=zeros(20,VesNum);
VesParam(1,:)=Len'*1e-6*40;        % Length
VesParam(2,:)=Diam'*1e-6*10;           % Diameter
VesParam(3,:)=WallTh'*1e-6;         % Wall thickness
% Young's modulus
for i=1:VesNum
  if VesCatogory(i)==1
    VesParam(4,i)=3.5e5;
  elseif VesCatogory(i)==3
    VesParam(4,i)=3.88e5;
  else
    VesParam(4,i)=3.7e5;
  end
end
VesParam(5,:)=1.26*ones(1,VesNum);  % Alpha
VesParam(6,:)=Visc'*1e-5;           % Viscosity
VesParam(7,:)=Hd;                   % Discharge hematocrit
VesParam(8,:)=0;                    % Venous pressure
VesParam(9,:)=3.25E-11;             % Peak of the input flow rate
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
VesParam(16,:)=10;                   % L order
VesParam(17,:)=10;                   % q order
if NoDim
  VesParam(18,:)=1e-3;                % scale_lamda
  VesParam(19,:)=1e-1;                % scale_u0
  VesParam(20,:)=1e-4;                % scale_r0
else
  VesParam(18,:)=1;                % scale_lamda
  VesParam(19,:)=1;                % scale_u0
  VesParam(20,:)=1;                % scale_r0
end

% inFileName=['Single_q' int2str(VesParam(16,1)) '_L' int2str(VesParam(15,1)) '_IN.bcs'];
inFileName=['Single_IN.bcs'];
outFileName=['Single_OUT.bcs'];
% Scale the input velocity profile
MeanBottomARatio=VelProc(dt,inFileName,outFileName,Vel(1),Vel(1),Nstep*dt/0.8,VesParam(18,1));
% VelVec=[Vel(1)*100 1 0];SinVelProc(dt,inFileName,outFileName,VelVec,Nstep*dt/0.8,1);

% Generate .in file of the simulation
fileName=GenInput(VesType,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);
system('copy Single*.* E:\Projects\Cpp\1DWin_MKL\1DWin_MKL\');
system('copy Single*.* E:\1DWin\Single\');

% system(['1DWin_MKL_IMP_NODIM.exe ' fileName '&']);
