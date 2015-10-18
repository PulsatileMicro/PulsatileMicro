%% Pries Hemodynamic Model
% Ref. Pries, A.R., et al., Blood flow in microvascular networks. Experiments and simulation[J], Circulation Research, 1990. 67(4): 826-834.
clear;clc;close all;
%% Ѫ���������� %%%%%%%% Ѫ���������� %%%%
% �궨��
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
Net_546_ID=1;       % 1 - 546��ϵĤѪ�����磨����Ӧ��: Net_546
Net_546_Meas_ID=2;  % 2 - 546��ϵĤѪ�����磨����Ӧǰ����Net_546_Meas
Egg_818_ID=3;       % 3 - 818����Ѫ������: Egg_818
Net_122_ID=4;       % 4 - 122�˹�Ѫ������: Net_122
Net_389_ID=5;       % 5 - 389��ϵĤ���磺Net_389
Net_913_ID=6;       % 6 - 913��ϵĤ���磺Net_913
Egg_CAM_ID=7;       % 7 - CAM����Ѫ�����磺Egg_CAM
Sub_CAM_ID=8;       % 8 - CAM����Ѫ�����������磺Sub_CAM
Egg_636_ID=9;       % 9 - 636����Ѫ�����磺Egg_636

VesType=1;
% ����΢ѭ����������Ϣ�͹�����Ϣ
switch VesType
  case Net_546_ID
    DatFile='T2810_h.dat';
    PrnFile='T2810_h_adapted.prn';
  case Net_546_Meas_ID
    DatFile='T2810_h.dat';
    PrnFile='';
  case Egg_818_ID
    DatFile='Egg818.dat';
    PrnFile='';
  case Net_122_ID
    DatFile='test.dat';
    PrnFile='test.prn';
  case Net_389_ID
    DatFile='T15_2_h.dat';
    PrnFile='';
  case Net_913_ID
    DatFile='1_8_TWS_morph.dat';
    PrnFile='';
  case Egg_CAM_ID
    % TODO:ԭʼ�͸��º���ļ���������ͬ
    DatFile='CAM_morph_update.dat';
    PrnFile='';
  case Sub_CAM_ID
    DatFile='SubCAM.dat';
    PrnFile='';
  case Egg_636_ID
    DatFile='571_Morph3_morph.dat';
    PrnFile='';
end
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);

%% Damping�����ų�ģʽ: DampFactor
% 1 - ��ʼģʽ
% 2 - ����ģ������10��
% 3 - ճ�Ͷȼ�С10��
% 4 - ����ģ����С10��
% 5 - ճ�Ͷ�����10��
DampFactor=1;

VesNum=length(DataArray(:,1));
SegName=DataArray(:,1);
SegType=DataArray(:,2);
From=DataArray(:,3);
To=DataArray(:,4);
Len=DataArray(:,6);
if VesType~=Net_546_ID && VesType~=Net_122_ID
  Diam=DataArray(:,5);
  Visc=2*ones(VesNum,1);
  if DampFactor==3
    Visc=Visc/10;
  elseif DampFactor==5
    Visc=Visc*10;
  end
else
  Diam=FuncPara(:,3);
  Visc=FuncPara(:,9);
  if DampFactor==3
    Visc=Visc/10;
  elseif DampFactor==5
    Visc=Visc*10;
  end
  WallTh=FuncPara(:,20);
end

% TempMod
% Visc(SegType==2)=Visc(SegType==2)/10;
% Visc=Visc/10;
% TempMod


% AdjMode �ܾ�����ģʽ
% 0: ���ֲ���
% 1�����ԷŴ�2��ָ���Ŵ�3�������Ŵ�
% 4��������С��5��ָ����С��6��������С
[Diam DiamRatio]=AdjustDiam(Diam,SegType,0);

% ����CAM������ָ�ѹ����ԭ��
if VesType==Egg_CAM_ID
  Boundary(10:11,3)=0.35*Boundary(10:11,3);
end
% �߽�����
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);
%%%%%%%
% �����߽�����
for i=1:length(BoundType)
  if BoundType(i)==1
    BoundFlow(i)=BoundFlow(i);
%     if BoundFlow(i)>0
%       BoundFlow(i)=BoundFlow(i).*DiamRatio(BoundNode(i)==From);
%     else
%       BoundFlow(i)=BoundFlow(i).*DiamRatio(BoundNode(i)==To);
%     end
  end
end
%%%%%%%

% �߽�Ѫ�ܷ������� %%%%
% ���ݱ߽������ķ�������
% Ѫ��Ϊ����->��߽�->�߽�ڵ���From����
% Ѫ��Ϊ����->���߽�->�߽�ڵ���To����
[From,To]=OriNodeModify(Boundary,From,To);

Qref=0.001; % TODO: Qrefû�ã����������Ҫ�ع�
% Boundary(:,4)=0.45;
BHd=BoundaryInput(Boundary,From,To,Qref);

% ��λ����
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s

% ѭ��1��ճ�Ͷȷ���
T=0;
if VesType==Net_546_ID || VesType==Net_122_ID
  % ���������Ӧ���ݣ���ôֻ��һ��
  CTime=1;
else
  % ���ѭ������ TODO(panqing):�������ж�
  CTime=5;
end
DebugVisc=zeros(VesNum,CTime);
DebugHd=zeros(VesNum,CTime);
DebugPressure=zeros(VesNum,CTime);
DebugFlow=zeros(VesNum,CTime);
Porder=1:VesNum;Norder=VesNum:1;  % ��ʼHd����˳��

while T<CTime   % ѭ��1��Visc����
  T=T+1
  % ���Է������ģ��
  [MeanP,DeltaP,MeanFlow]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam,Visc,Len,From,To,DiamRatio);
  % ��������ģ��
  [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP);
  %   FromNew=From;ToNew=To;MeanFlowNew=MeanFlow;DeltaPNew=DeltaP;
  % ����Hd����˳��
  [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,2,0);
  % ����Hd
  [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam,MeanFlowNew,0);
  
  if VesType~=Net_546_ID && VesType~=Net_122_ID
    %%ճ�Ͷȷ���
    umDiam=Diam.*1e6;   %um
    for i=1:VesNum
      Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
      if DampFactor==3
        Visc(i)=Visc(i)/10;
      elseif DampFactor==5
        Visc(i)=Visc(i)*10;
      end
    end
    if T>1
      Visc=0.5*(Visc+DebugVisc(:,T-1));
    end
  end
  
  % ��¼ÿ�ε����ķ�����
  DebugVisc(:,T)=Visc;
  DebugHd(:,T)=Hd;
  DebugPressure(:,T)=MeanP;
  DebugFlow(:,T)=MeanFlowNew;
end
Vel=MeanFlowNew/1e12/60./(pi*0.25*Diam.^2)*1e3;

if VesType==Net_913_ID
  SegType=AutoVesType(Boundary,FromNew,ToNew,Porder);
elseif VesType==Egg_818_ID
  ConstSegType=zeros(VesNum,1);
  % ����Ѫ�ܶ��޷�ͨ���Զ��㷨�ж�����
  ConstSegType([216;544;547;559;561;47;57;70;80;143;681;683])=[3;ones(11,1)];
  SegType=autoVesTypeAdv(ConstSegType,FromNew,ToNew,Porder);
elseif VesType==Egg_636_ID
  ConstSegType=zeros(VesNum,1);
  ConstSegType([167;168;143;144;560;561;562;212;580;581;583;410;615;449;514;452])=ones(16,1);
  ConstSegType([213;544;506;507;583;394;384;385;386;396])=2*ones(10,1);
  SegType=autoVesTypeAdv(ConstSegType,FromNew,ToNew,Porder);
  SegType(424)=3;
  SegType(428)=2;
end

%% ����Ƚ�
% Ѫѹ��Ѫ��
figure;
if VesType==Net_546_ID
  subplot(511);
  plot(1:VesNum,FuncPara(:,8),'-ro',1:VesNum,MeanFlow,'-b.');title('Flow Rate');
  subplot(512);
  plot(1:VesNum,FuncPara(:,4),'-ro',1:VesNum,MeanP,'-b.');title('Pressure');
  subplot(513)
  plot(1:VesNum,FuncPara(:,9),'-ro',1:VesNum,Visc*1e3,'-b.');title('Viscosity');
  subplot(514)
  plot(1:VesNum,FuncPara(:,24),'-ro',1:VesNum,DeltaP,'-b.');title('DeltaP');
  subplot(515)
  plot(1:VesNum,FuncPara(:,7),'-ro',1:VesNum,Vel,'-b.');title('Vel');
else
  subplot(211);plot(MeanFlowNew);title('Flow Rate');
  subplot(212);plot(MeanP);title('Pressure');
end

% ���Ʋ�������ͼ
PlotTopol(MeanP,VesType);
% PlotTopol(MeanFlowNew,VesType);
% PlotTopol(SegType,VesType);

%% �������
Visc=Visc*1e3;
Flow=MeanFlowNew;
DeltaP=DeltaPNew;
% ����steady state����
switch VesType
  case Net_546_ID
    if DampFactor==3
      save('546_Adap_Visc-1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    elseif DampFactor==5
      save('546_Adap_Visc+1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    else
      save('546_Adap.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    end
  case Net_546_Meas_ID
    if DampFactor==3
      save('546_Meas_Visc-1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    elseif DampFactor==5
      save('546_Meas_Visc+1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    else
      save('546_Meas.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    end
  case Egg_818_ID
    if DampFactor==3
      save('Egg818_Meas_Visc-1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    elseif DampFactor==5
      save('Egg818_Meas_Visc+1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    else
      save('Egg818_Meas.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    end
  case Net_389_ID
    if DampFactor==3
      save('389_Meas_Visc-1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    elseif DampFactor==5
      save('389_Meas_Visc+1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    else
      save('389_Meas.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    end
  case Net_913_ID
    if DampFactor==3
      save('913_Meas_Visc-1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    elseif DampFactor==5
      save('913_Meas_Visc+1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    else
      save('913_Meas.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    end
  case Sub_CAM_ID
    if DampFactor==3
      save('SubCAM_Meas_Visc-1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    elseif DampFactor==5
      save('SubCAM_Meas_Visc+1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    else
      save('SubCAM_Meas.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    end
  case Egg_636_ID
    if DampFactor==3
      save('Egg_636_Meas_Visc-1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    elseif DampFactor==5
      save('Egg_636_Meas_Visc+1.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    else
      save('Egg_636_Meas.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
    end
end

%% ���±���.dat�ļ�
% 1Dģ�ͳ���Ҫ��Ѫ������������ܳ�������������������
% ���·���Ѫ������
switch VesType
  case Net_546_ID
    fiddat=fopen('Net_546_Adap.dat','w');
    % Header
    fprintf(fiddat,'RAT MESENTERY  01-08-85\n');
    fprintf(fiddat,'546\n');
    fprintf(fiddat,'\n\n\n');
    % Content
    fprintf(fiddat,'SEGMENT NAME   TYPE   FROM  TO    DIAM   LENGTH  HDmes  VELmes\n');
    for i=1:length(SegName)
      fprintf(fiddat,'%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n',SegName(i),SegType(i),FromNew(i),ToNew(i),Diam(i)*1e6,Len(i)*1e6,Hd(i),0);
    end
    % Boundary
    fprintf(fiddat,'%d\n',36);
    fprintf(fiddat,'INOD  IPRFL  PR/FL  HEMAT  IPRFL=0 FOR PRESSURE, 1 FOR FLOW B.C.\n');
    for i=1:36
      fprintf(fiddat,'%d\t%d\t%f\t%f\n',BoundNode(i),BoundType(i),BoundFlow(i),BoundHd(i));
    end
    fclose(fiddat);
    fidprn=fopen('Net_546_Adap.prn','w');
    fprintf(fidprn,'Name  oldDia   newDia   Pr_mean  Tau_real  Tau_fit    Vel       Flow    Visc      HD     P_O2      S_O2     S_tot   MetStim   Hydrod   CondSti MetLoc   MetConv   Typ   WallTh   Sigma    SO2in    SO2out     deltaP   ConsAdFact\n\n\n\n');
    for i=1:length(SegName)
      fprintf(fidprn,'%d\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',...
        SegName(i),Diam(i)*1e6,Diam(i)*1e6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,WallTh(i).*DiamRatio(i),0,0,0,0,0);
    end
    fclose(fidprn);
  case Net_546_Meas_ID
    fiddat=fopen('T2810_meas.dat','w');
    % Header
    fprintf(fiddat,'RAT MESENTERY  01-08-85\n');
    fprintf(fiddat,'546\n');
    fprintf(fiddat,'\n\n\n');
    % Content
    fprintf(fiddat,'SEGMENT NAME   TYPE   FROM  TO    DIAM   LENGTH  HDmes  VELmes\n');
    for i=1:length(SegName)
      fprintf(fiddat,'%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n',SegName(i),SegType(i),FromNew(i),ToNew(i),Diam(i)*1e6,Len(i)*1e6,Hd(i),0);
    end
    % Boundary
    fprintf(fiddat,'%d\n',36);
    fprintf(fiddat,'INOD  IPRFL  PR/FL  HEMAT  IPRFL=0 FOR PRESSURE, 1 FOR FLOW B.C.\n');
    for i=1:36
      fprintf(fiddat,'%d\t%d\t%f\t%f\n',BoundNode(i),BoundType(i),BoundFlow(i),BoundHd(i));
    end
    fclose(fiddat);
  case Egg_818_ID
    fiddat=fopen('Egg818_update.dat','w');
    fprintf(fiddat,'Chicken Embryo CAM\n');
    fprintf(fiddat,'818 XXX NSEG NNOD\n');
    fprintf(fiddat,'\n\n\n'); % ��֮ǰ��ʽ����һ��
    fprintf(fiddat,'SEGMENT NAME   TYPE   FROM       TO    DIAM   LENGTH  HDmes  VELmes\n');
    for i=1:length(SegName)
      fprintf(fiddat,'%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n', SegName(i),SegType(i),FromNew(i),ToNew(i),Diam(i)*1e6,Len(i)*1e6,Hd(i),0);
    end
    % Boundary
    fprintf(fiddat,'%d\n',43);
    fprintf(fiddat,'INOD  IPRFL  PR/FL  HEMAT  IPRFL=0 FOR PRESSURE, 1 FOR FLOW B.C.\n');
    for i=1:43
      fprintf(fiddat,'%d\t%d\t%f\t%f\n',BoundNode(i),BoundType(i),BoundFlow(i),BoundHd(i));
    end
    fclose(fiddat);
  case Net_389_ID
    fiddat=fopen('T15_2_update.dat','w');
    % Header
    fprintf(fiddat,'RAT MESENTERY  01-08-85\n');
    fprintf(fiddat,'389\n');
    fprintf(fiddat,'\n\n\n');
    % Content
    fprintf(fiddat,'SEGMENT NAME   TYPE   FROM  TO    DIAM   LENGTH  HDmes  VELmes\n');
    for i=1:length(SegName)
      fprintf(fiddat,'%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n',SegName(i),SegType(i),FromNew(i),ToNew(i),Diam(i)*1e6,Len(i)*1e6,Hd(i),0);
    end
    % Boundary
    fprintf(fiddat,'%d\n',22);
    fprintf(fiddat,'INOD  IPRFL  PR/FL  HEMAT  IPRFL=0 FOR PRESSURE, 1 FOR FLOW B.C.\n');
    for i=1:22
      fprintf(fiddat,'%d\t%d\t%f\t%f\n',BoundNode(i),BoundType(i),BoundFlow(i),BoundHd(i));
    end
    fclose(fiddat);
  case Net_913_ID
    fiddat=fopen('1_8_TWS_morph_update.dat','w');
    % Header
    fprintf(fiddat,'RAT MESENTERY  01-08-85\n');
    fprintf(fiddat,'913\n');
    fprintf(fiddat,'\n\n\n');
    % Content
    fprintf(fiddat,'SEGMENT NAME   TYPE   FROM  TO    DIAM   LENGTH  HDmes  VELmes\n');
    for i=1:length(SegName)
      fprintf(fiddat,'%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n',SegName(i),SegType(i),FromNew(i),ToNew(i),Diam(i)*1e6,Len(i)*1e6,Hd(i),0);
    end
    % Boundary
    fprintf(fiddat,'%d\n',65);
    fprintf(fiddat,'INOD  IPRFL  PR/FL  HEMAT  IPRFL=0 FOR PRESSURE, 1 FOR FLOW B.C.\n');
    for i=1:65
      fprintf(fiddat,'%d\t%d\t%f\t%f\n',BoundNode(i),BoundType(i),BoundFlow(i),BoundHd(i));
    end
    fclose(fiddat);
  case Egg_CAM_ID
    fiddat=fopen('CAM_morph_update.dat','w');
    fprintf(fiddat,'Chicken Embryo CAM\n');
    fprintf(fiddat,'7128 XXX NSEG NNOD\n');
    fprintf(fiddat,'\n\n\n'); % ��֮ǰ��ʽ����һ��
    fprintf(fiddat,'SEGMENT NAME   TYPE   FROM       TO    DIAM   LENGTH  HDmes  VELmes\n');
    for i=1:length(SegName)
      fprintf(fiddat,'%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n', SegName(i),SegType(i),FromNew(i),ToNew(i),Diam(i)*1e6,Len(i)*1e6,Hd(i),0);
    end
    % Boundary
    fprintf(fiddat,'%d\n',48);
    fprintf(fiddat,'INOD  IPRFL  PR/FL  HEMAT  IPRFL=0 FOR PRESSURE, 1 FOR FLOW B.C.\n');
    for i=1:48
      fprintf(fiddat,'%d\t%d\t%f\t%f\n',BoundNode(i),BoundType(i),BoundFlow(i),BoundHd(i));
    end
    fclose(fiddat);
  case Egg_636_ID
    fiddat=fopen('Egg636_morph_update.dat','w');
    fprintf(fiddat,'Chicken Embryo CAM\n');
    fprintf(fiddat,'636 XXX NSEG NNOD\n');
    fprintf(fiddat,'\n\n\n'); % ��֮ǰ��ʽ����һ��
    fprintf(fiddat,'SEGMENT NAME   TYPE   FROM       TO    DIAM   LENGTH  HDmes  VELmes\n');
    for i=1:length(SegName)
      fprintf(fiddat,'%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n', SegName(i),SegType(i),FromNew(i),ToNew(i),Diam(i)*1e6,Len(i)*1e6,Hd(i),0);
    end
    % Boundary
    fprintf(fiddat,'%d\n',45);
    fprintf(fiddat,'INOD  IPRFL  PR/FL  HEMAT  IPRFL=0 FOR PRESSURE, 1 FOR FLOW B.C.\n');
    for i=1:45
      fprintf(fiddat,'%d\t%d\t%f\t%f\n',BoundNode(i),BoundType(i),BoundFlow(i),BoundHd(i));
    end
    fclose(fiddat);
end
