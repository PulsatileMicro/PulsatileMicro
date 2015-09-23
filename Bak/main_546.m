%% Main function
clear;clc;close all;

VesType='Net_546';        % Network name
DatFile='T2810_h.DAT';    % Network data file name

ModelParam=zeros(12,1);
% Measured/Adapted switcher 0: use adapted diameter, 1: use measured diameter
ModelParam(1)=0;
MeasDiam=ModelParam(1);
% Nondimensionalization flag 1: Nondimensionalization
ModelParam(2)=0;
NoDim=ModelParam(2);
% Viscosity update strategy 0: No viscosity update 1: vitro 2: ESL vivo 3: No ESL vivo
ModelParam(3)=2;
% Riemann solver flag 1: Use riemann solver for the BioFlux
ModelParam(4)=0;
% Binary output flag 0: Output text format 1: Output binary format
ModelParam(5)=0;
% Time step
ModelParam(6)=1e-3;
dt=ModelParam(6);
% Number of steps
Period=0.8;
% Period=0.325;
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
if MeasDiam
  PrnFile='T2810_h_ini.prn';
else
  PrnFile='T2810_h_adapted.prn';
end
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);

SegName=DataArray(:,1);
From=DataArray(:,3);
To=DataArray(:,4);
Len=DataArray(:,6);
VesCategory=DataArray(:,2);
VesNum=length(SegName);

if MeasDiam
  % Measured data
  Diam=FuncPara(:,3);    % Measured
  Hd=FuncPara(:,10);    % Predicted
  %     DeltaP=FuncPara(:,24);
  WallTh=Diam*0.1;
  load 546_Meas.mat;
else
  % Predicted data
  Diam=FuncPara(:,3);   % Adapted
  WallTh=FuncPara(:,20);
  %     Vel=FuncPara(:,7);    % Predicted
  load 546_Adap.mat;
end

BCTypeAll=[];
Bifur=zeros(VesNum,4);
BCVal=zeros(VesNum,1);
BoundVel=[];
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
    if From(i) == 830
%       BCType=[BCType 'u'];
      BCType=[BCType 'q'];
      BCVal(i)=Vel(i);
%       BCType=[BCType 'p'];
%       BCVal(i)=MeanP(i);
      Bifur(i,1:2)=[0.4338 2];
      BoundVel=[BoundVel;BCVal(i) i 0];
    else
      % pulsatile u input
%       BCType=[BCType 'u'];
      BCType=[BCType 'q'];
%       BCType=[BCType 'p'];
      tmpInd=find(Boundary(:,1)==From(i));
      Bifur(i,1:2)=[Boundary(tmpInd,4) 2];
      if MeasDiam
%         BCVal(i)=Boundary(tmpInd,3)/60/1e12/(pi*(Diam(i)*1e-6).^2/4);
        BCVal(i)=Vel(i);
%         BCVal(i)=MeanP(i);
      else
        BCVal(i)=Vel(i);
%         BCVal(i)=MeanP(i);
      end
      BoundVel=[BoundVel;BCVal(i) i 0];
    end
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
    if To(i)==825
%       BCType=[BCType 'T'];
%       BCVal(i)=0;
      BCType=[BCType 'R'];
      BCVal(i)=(MeanP(i)-DeltaP(i)/2)/Flow(i)*133*60*1e12;
    else
      BCType=[BCType 'R'];
      BCVal(i)=(MeanP(i)-DeltaP(i)/2)/Flow(i)*133*60*1e12;
    end
    Bifur(i,3:4)=[0 0];
    % BCType=[BCType 'T'];
    % BCVal(i)=0;
    
    % u boundary
    %         BCType=[BCType 'u'];
    %         tmpInd=find(Boundary(:,1)==To(i));
    %         Bifur(i,3:4)=[Boundary(tmpInd,4) 2];
    %         BCVal(i)=Vel(i);
    %         BoundVel=[BoundVel;BCVal(i) i 1];
  else
    BCType=[BCType 'J'];
    Bifur(i,3:4)=outInd1;
  end
  BCTypeAll=[BCTypeAll;BCType];
end

VesParam=zeros(20,VesNum);
VesParam(1,:)=Len'*1e-6;            % Length
VesParam(2,:)=Diam'*1e-6;           % Diameter
VesParam(3,:)=WallTh'*1e-6;         % Wall thickness
for i=1:VesNum
  if VesCategory(i)==1
    VesParam(4,i)=3.5e5;
  elseif VesCategory(i)==3
    VesParam(4,i)=3.88e5;
%     VesParam(4,i)=5.4e6;
  else
    VesParam(4,i)=3.7e5;
%     VesParam(4,i)=1.8e6;
  end
end
VesParam(5,:)=1.26*ones(1,VesNum);  % Alpha
VesParam(6,:)=Visc'*1e-3;           % Viscosity
% VesParam(6,:)=Visc'*1e-4;           % Viscosity
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
VesParam(16,:)=3;                   % q
VesParam(17,:)=3;                   % L
if NoDim
  VesParam(18,:)=1e-3;                % scale_lamda
  VesParam(19,:)=1e-1;                % scale_u0
  VesParam(20,:)=1e-4;                % scale_r0
else
  VesParam(18,:)=1;
  VesParam(19,:)=1;
  VesParam(20,:)=1;
end
VesParam(21,:)=3;                   % Number of history points
VesParam(22,:)=0;                   % Taper Rate

% Scale the input velocity profile
inFileName='Net_546_IN_';
outFileName='Net_546_OUT_';
% MeanBottomRatioA=VelProc546_Gaehtgens(dt,inFileName,outFileName,BoundVel,Nstep*dt/Period,VesParam(18,1));
% MeanBottomRatioA=FlowProc546_Gaehtgens(dt,inFileName,outFileName,BoundVel,Nstep*dt/Period,Diam,VesParam(18,1),VesParam(19,1));
% MeanBottomRatioA=VelProc546(dt,inFileName,outFileName,BoundVel,Nstep*dt/Period,VesParam(18,1));
MeanBottomRatioA=FlowProc546(dt,inFileName,outFileName,BoundVel,Nstep*dt/Period,Diam,VesParam(18,1),VesParam(19,1));
% MeanBottomRatioA=PressureProc546(dt,inFileName,outFileName,BoundVel,Nstep*dt/Period,Diam,VesParam(18,1),VesParam(19,1));
% SinVelProc(dt,inFileName,outFileName,BoundVel,Nstep*dt/Period,VesParam(18,1));MeanBottomRatioA=1;
ModelParam(8)=MeanBottomRatioA;

% Generate .in file of the simulation
fileName=GenInput(VesType,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);
system('copy Net_546*.* E:\Projects\CPP\1DWin_MKL\1DWin_MKL\');
system('copy Net_546*.* E:\1DWin\546\');