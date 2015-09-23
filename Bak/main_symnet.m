clear;clc;close all;
order=5;    % order指从起点到毛细血管层的阶数，总阶数为2*order-1
curOrder=0;
cnt=0;
lastNum=0;
VesNum=2^order-1 + 2^(order-1)-1;

% Parameters for PI reduction 60%
% startDiam=4e-5
% ratio_L=50
% ratio_D_L=0.85
% Visc=5e-3
% T=0.8
% E=1e5

% Parameters for PI reduction 30%
% startDiam=4e-5
% ratio_L=60;
% ratio_D_L=0.85
% Visc=5e-3
% T=0.8
% E=3.5e5

% Parameters of start vessel
% startDiam=4e-5/0.68;
startDiam=4e-5;
ratio_L=70;
ratio_D_L=0.89;
ratio_D_R=0.89;
ratio_D_Conv=0.89;
% ratio_D_L=0.7583;
% ratio_D_R=0.7583;
% ratio_D_Conv=0.7583;

DiamVec=zeros(1,VesNum);
LenVec=zeros(1,VesNum);
VesNumVec=zeros(1,VesNum);
Bifur=zeros(4,VesNum);
BCTypeAll=[];
SegName=zeros(1,VesNum);

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
      BCType=[BCType 'B' 'C'];
    elseif i == 1
      BCType=[BCType 'u' 'B'];
%       BCType=[BCType 'p' 'B'];
    else
      BCType=[BCType 'B' 'B'];
    end
  else
    if cnt==2^(2*order-curOrder-1)
      cnt=0;
      curOrder=curOrder+1;
    end
    lastNum=2*2^(2*order-curOrder-1);
    DiamVec(i)=DiamVec(i-lastNum+cnt)/ratio_D_Conv;
    LenVec(i)=DiamVec(i)*ratio_L;
    Bifur(1,i)=i-lastNum+cnt;
    Bifur(2,i)=i-lastNum+cnt+1;
    Bifur(3,i-lastNum+cnt)=i;
    Bifur(4,i-lastNum+cnt)=i-lastNum+cnt+1;
    Bifur(3,i-lastNum+cnt+1)=i-lastNum+cnt;
    Bifur(4,i-lastNum+cnt+1)=i;
    cnt=cnt+1;
    if i==VesNum
      BCType=[BCType 'C' 'T'];
      BCVal(i)=0.8;
    else
      BCType=[BCType 'C' 'C'];
    end
  end
  BCTypeAll=[BCTypeAll;BCType];
end

LenVec=[0.00280000000000000,0.00249200000000000,0.00249200000000000,0.00221788000000000,0.00221788000000000,0.00221788000000000,0.00221788000000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00175678274800000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00197391320000000,0.00221788000000000,0.00221788000000000,0.00221788000000000,0.00221788000000000,0.00249200000000000,0.00249200000000000,0.00280000000000000;];
Vel=1;
VesParam=zeros(22,VesNum);
VesParam(1,:)=LenVec;
VesParam(2,:)=DiamVec;
VesParam(3,:)=0.1*startDiam;
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
VesParam(21,:)=3;                   % Number of history points
VesParam(22,:)=0;                   % Taper Rate

VesType='SymNet';
% Time step
dt=1e-3;
% Number of step
Nstep=4.8/dt;
inFileName='SymNet_IN_1.bcs';
outFileName='SymNet_OUT_1.bcs';
% Scale the input velocity profile
VelProc(dt,inFileName,outFileName,Vel,Vel,Nstep*dt/0.8,1);
% PressureProc(dt,inFileName,outFileName,100,100,Nstep*dt/0.8,1);
% VelVec=[Vel(1) 1 0];SinVelProc(dt,inFileName,outFileName,VelVec,Nstep*dt/0.8,1);

fileName=GenSymNet(VesType, VesParam, BCTypeAll, Bifur', dt, Nstep);
system('copy SymNet*.* E:\Projects\CPP\1DWin_MKL\1DWin_MKL\');
system('copy SymNet*.* E:\1DWin\SymNet\');