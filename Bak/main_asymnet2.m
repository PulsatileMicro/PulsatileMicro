% Asymmetric network type II
% ----------
%    |  |  |
% ----------
clear;clc;close all;
order=2;
VesNum=1+3*order;

% Parameters of start vessel
startDiam=4e-3;
ratio_L=120;
% ratio_L=20;
% ratio_L=50*0.85/0.77*0.74/0.70;
ratio_D_L=0.85;
ratio_D_R=0.85;
% ratio_D_L=0.75;
% ratio_D_R=0.75;

DiamVec=zeros(1,VesNum);
LenVec=zeros(1,VesNum);
VesNumVec=zeros(1,VesNum);
Bifur=zeros(4,VesNum);
BCTypeAll=[];
SegName=zeros(1,VesNum);

Bifur=[0 2 2 3;1 3 6 7;1 2 4 5;3 5 5 6;3 4 4 6;4 5 2 7;2 6 0 0]';
BCTypeAll=['u' 'B';'B' 'C';'B' 'B';'B' 'C';'B' 'C';'C' 'C';'C' 'T'];
DiamVec(1)=startDiam;
DiamVec(2)=DiamVec(1)*ratio_D_L^2;
DiamVec(3)=DiamVec(1)*ratio_D_L;
DiamVec(4:5)=DiamVec(3)*ratio_D_L;
DiamVec(6)=DiamVec(4)/ratio_D_L;
DiamVec(7)=DiamVec(2)/ratio_D_L^2;
LenVec=DiamVec*ratio_L;

% Vel=5e4;
Vel=10;
VesParam=zeros(15,VesNum);
VesParam(1,:)=LenVec;
VesParam(2,:)=DiamVec;
VesParam(3,:)=0.4*DiamVec/2;         % Wall thickness
VesParam(4,:)=3.5e5*ones(1,VesNum);  % Young's modulus
VesParam(5,:)=1.26*ones(1,VesNum);  % Alpha
VesParam(6,:)=5e-3;                 % Viscosity
VesParam(7,:)=0.45;                 % Discharge hematocrit
VesParam(8,:)=0;                    % Venous pressure
VesParam(9,:)=3.25E-11;             % Peak of the input flow rate
VesParam(10,:)=0;                   % Segment name
VesParam(11,:)=0;                   % From nodes
VesParam(12,:)=0;                   % To nodes
VesParam(13,:)=Vel;

VesType='ASymNet2';
% Time step
dt=1e-3;
% Number of step
Nstep=4.8/dt;
inFileName='ASymNet2_IN_1.bcs';
outFileName='ASymNet2_OUT_1.bcs';
% Scale the input velocity profile
VelProc(dt,inFileName,outFileName,Vel,Vel,Nstep*dt/0.8,1);
% PressureProc(dt,inFileName,outFileName,100,10,Nstep*dt/0.8,1);
% VelVec=[Vel(1) 1 0];SinVelProc(dt,inFileName,outFileName,VelVec,Nstep*dt/0.8,1);

% VesParam(14,:)=Eval_Gamma(VesParam(2,:),1e4*ones(1,VesNum),VesParam(3,:));  % Visc part of viscoelasticity
% VesParam(15,:)=Eval_GammaII(VesParam(2,:),1e6*ones(1,VesNum));
fileName=GenSymNet(VesType, VesParam, BCTypeAll, Bifur', dt, Nstep);
system('copy ASymNet2*.* E:\Projects\CPP\1DWin_MKL\1DWin_MKL\');
system('copy ASymNet2*.* E:\1DWin\SymNet\');