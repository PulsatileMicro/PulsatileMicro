%% Pries Hemodynamic Model
% Ref. Pries, A.R., et al., Blood flow in microvascular networks. Experiments and simulation[J], Circulation Research, 1990. 67(4): 826-834.
clear;clc;close all;
%% 血管网络类型 %%%%%%%% 血管网络类型 %%%%
% 宏定义
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
Net_546_ID=1;       % 1 - 546肠系膜血管网络（自适应后）: Net_546
Net_546_Meas_ID=2;  % 2 - 546肠系膜血管网络（自适应前）：Net_546_Meas
Egg_818_ID=3;       % 3 - 818鸡胚血管网络: Egg_818
Net_122_ID=4;       % 4 - 122人工血管网络: Net_122
Net_389_ID=5;       % 5 - 389肠系膜网络：Net_389
Net_913_ID=6;       % 6 - 913肠系膜网络：Net_913
Egg_CAM_ID=7;       % 7 - CAM鸡胚血管网络：Egg_CAM
Sub_CAM_ID=8;       % 8 - CAM鸡胚血管网络子网络：Sub_CAM
Egg_636_ID=9;       % 9 - 636鸡胚血管网络：Egg_636

VesType=1;
% 导入微循环的拓扑信息和功能信息
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
    % TODO:原始和更新后的文件名不能相同
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

%% Damping因素排除模式: DampFactor
% 1 - 初始模式
% 2 - 杨氏模量增大10倍
% 3 - 粘滞度减小10倍
% 4 - 杨氏模量减小10倍
% 5 - 粘滞度增大10倍
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


% AdjMode 管径调节模式
% 0: 保持不变
% 1：线性放大，2：指数放大，3：对数放大
% 4：线性缩小，5：指数缩小，6：对数缩小
[Diam DiamRatio]=AdjustDiam(Diam,SegType,0);

% 查找CAM网络出现负压降的原因
if VesType==Egg_CAM_ID
  Boundary(10:11,3)=0.35*Boundary(10:11,3);
end
% 边界数据
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);
%%%%%%%
% 调整边界流量
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

% 边界血管方向修正 %%%%
% 依据边界条件的符号修正
% 血流为正数->入边界->边界节点在From类中
% 血流为负数->出边界->边界节点在To类中
[From,To]=OriNodeModify(Boundary,From,To);

Qref=0.001; % TODO: Qref没用，这个程序需要重构
% Boundary(:,4)=0.45;
BHd=BoundaryInput(Boundary,From,To,Qref);

% 单位调整
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s

% 循环1，粘滞度反馈
T=0;
if VesType==Net_546_ID || VesType==Net_122_ID
  % 如果有自适应数据，那么只跑一次
  CTime=1;
else
  % 最大循环次数 TODO(panqing):以收敛判断
  CTime=5;
end
DebugVisc=zeros(VesNum,CTime);
DebugHd=zeros(VesNum,CTime);
DebugPressure=zeros(VesNum,CTime);
DebugFlow=zeros(VesNum,CTime);
Porder=1:VesNum;Norder=VesNum:1;  % 初始Hd计算顺序

while T<CTime   % 循环1，Visc反馈
  T=T+1
  % 线性方程求解模块
  [MeanP,DeltaP,MeanFlow]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam,Visc,Len,From,To,DiamRatio);
  % 逆流处理模块
  [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP);
  %   FromNew=From;ToNew=To;MeanFlowNew=MeanFlow;DeltaPNew=DeltaP;
  % 调整Hd计算顺序
  [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,2,0);
  % 计算Hd
  [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam,MeanFlowNew,0);
  
  if VesType~=Net_546_ID && VesType~=Net_122_ID
    %%粘滞度反馈
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
  
  % 记录每次迭代的仿真结果
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
  % 部分血管段无法通过自动算法判断类型
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

%% 结果比较
% 血压、血流
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

% 绘制参数拓扑图
PlotTopol(MeanP,VesType);
% PlotTopol(MeanFlowNew,VesType);
% PlotTopol(SegType,VesType);

%% 结果保存
Visc=Visc*1e3;
Flow=MeanFlowNew;
DeltaP=DeltaPNew;
% 保存steady state数据
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

%% 重新保存.dat文件
% 1D模型程序要求血管连接情况不能出现零进三出或三进零出
% 重新分析血管类型
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
    fprintf(fiddat,'\n\n\n'); % 与之前格式保持一致
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
    fprintf(fiddat,'\n\n\n'); % 与之前格式保持一致
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
    fprintf(fiddat,'\n\n\n'); % 与之前格式保持一致
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
