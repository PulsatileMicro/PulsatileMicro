%% Generate input file of single artery for 1D program
clear;clc;close all;

% Network name
VesType='Artery_Loop';
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

PercentL=[0.8 0.9 1.0 1.1 1.2];
PercentD=[0.8 0.9 1.0 1.1 1.2];
PercentE=[0.8 0.9 1.0 1.1 1.2];
PercentU=[0.8 0.9 1.0 1.1 1.2];
PIDamp=zeros(5,1);
dA=zeros(5,1);

for i=1:length(PercentL)
  for k=1:length(PercentL)
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
    VesParam(2,:)=0.02*0.05*PercentD(i);        % Diameter
    VesParam(3,:)=0.001*0.05;         % Wall thickness
    VesParam(4,:)=0.6e6;          % Young's modulus
    VesParam(5,:)=1.1;            % Alpha
    VesParam(6,:)=4e-3;           % Viscosity
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
    VesParam(16,:)=10;            % q
    VesParam(17,:)=10;            % L
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
    
    Vel=1*PercentU(k);
    %   Vel=1;
    BCTypeAll=[];
    BCType=[];
    BCType=['u' 'T'];
    Bifur=[0 2];
    BCVal=0.0;
    % BCType=['p' 'T'];
    % Bifur=[0 2];
    % BCVal=0.0;
    BCTypeAll=[BCTypeAll;BCType];
    
    % Time step
    dt=1e-3;
    % Number of step
    Nstep=4.8/dt;
    inFileName='Artery_Loop_IN.bcs';
    outFileName='Artery_Loop_OUT.bcs';
    % Scale the input velocity profile
    % PressureProc(dt,inFileName,outFileName,150.8,0,Nstep*dt/0.8,1);   % Input Unit mmHg
    VelProc(dt,inFileName,outFileName,Vel,0,Nstep*dt/0.8,1);
    % FlowProc(dt,inFileName,outFileName,150.8/10,0,Nstep*dt/0.8,VesParam(2,:)*1e6,1,1);
    % VelVec=[2e4 1 0]; SinVelProc(dt,inFileName,outFileName,VelVec,Nstep*dt/0.8,1);
    
    % Generate .in file of the simulation
    fileName=GenInput(VesType,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);
    % Run the simulation
    winopen('run_artery_loop.bat');
    
    pause(3.5); % —” ±
    
    % Analyze the results
    PeriodNum=6;    % Number of cardiac cycle in the simulation
    t_plot=1+800*5:800*PeriodNum;
    all_plot=1:800*PeriodNum;
    NumHisPt=VesParam(21,1);
    [MeanP MeanQ MeanU MeanA MeanVisc tAll PAll UAll QAll AAll ViscAll]=HemoAnal4Network(NumHisPt,VesType,VesNum,t_plot,all_plot);
    R=8*MeanVisc(1,2).*VesParam(1)/pi./(((MeanA(1,2)+MeanA(1,2))/2/pi).^2);
    MeanP=MeanP/133;
    MeanQ=MeanQ*1e12*60;
    MeanU=MeanU*1000;
    MeanA=MeanA*1e12;
    PAll=PAll/133;
    QAll=QAll*1e12*60;
    UAll=UAll*1000;
    AAll=AAll*1e12;
    
    dA(i,k)=(MeanA(1)-MeanA(end))/mean(MeanA);
    for j=1:NumHisPt
      PI_U(j)=(max(UAll(j,t_plot))-min(UAll(j,t_plot)))./mean(UAll(j,t_plot));
      PI_P(j)=(max(PAll(j,t_plot))-min(PAll(j,t_plot)))./mean(PAll(j,t_plot));
    end
    PIDamp(i,k)=(PI_P(1)-PI_P(end))/PI_P(1);
    system('del Artery_Loop.his');
  end
end

