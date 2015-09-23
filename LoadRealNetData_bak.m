%% 根据仿真的血管网络类型，导入数据
%%%% 如果需要更改导入的数据，需修改这个函数 %%%%
function [ModelParam,DataArray,Boundary,FuncPara,DatMatrix,VesNum,TaperRate,PulsatLevel,DampFactorName]...
=LoadRealNetData(NetTypeID,ModelParam,DampFactor,RunType,ViscRatio,ERatio)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN E5 E02 Visc5 Visc02 DampE DampVisc
global ADAP NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE N_PERIOD DT PERIOD ODEREL ODEABS IN_BOUND_TYPE RANDOM_PHASE LIN_PA VEL_PROFILE BG_LEN_RATIO BG_MASS_RATIO USE_OPTBOUND
global ONED_PREP FREQ_ANAL BOUND_OPT_SS BOUND_OPT_PULSATILE

switch NetTypeID
    %% 546自适应后网络
  case Net_546_ID
    DatFile='Men_546.dat';
    PrnFile='Men_546.prn';
    %% 546自适应前网络
  case Net_546_Meas_ID
    DatFile='Men_546.dat';
    PrnFile='';
    %% 389血管网络
  case Net_389_ID
    DatFile='Men_389.dat';    % Network data file name
    PrnFile='';
    %% 913血管网络
  case Net_913_ID
    DatFile='Men_913.dat';    % Network data file name
    PrnFile='';
    %% 鸡胚818血管网络
  case Egg_818_ID
    DatFile='Vit_818.dat';
    PrnFile='';
    %% 鸡胚636血管网络
  case Egg_636_ID
    DatFile='Vit_636.dat';
    PrnFile='';
    %% CAM血管网络
  case Egg_CAM_ID
    DatFile='CAM_7128.DAT';    % Network data file name
    PrnFile='';
    %% SubCAM血管网络
  case Sub_CAM_ID
    DatFile='SubCAM.DAT';    % Network data file name
    PrnFile='';
    %% 122测试血管网络
  case Net_122_ID
    DatFile='Men_122.dat';
    PrnFile='Men_122.prn';
end

%% 对于所有真实血管网络，从文件读取数据
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
% 对于CAM网络，调节其边界流量值
% 因为边界流量值不是实测的，而是依据流速-管径的经验关系给出的，需要进行微调
% if NetTypeID==Egg_CAM_ID
%   Boundary(10:11,3)=0.35*Boundary(10:11,3);
% end

% 边界数据
if ModelParam(USE_OPTBOUND)==1
  switch NetTypeID
    case Net_546_Meas_ID
      load Men_546_Meas_BoundOpt.mat;
    case Net_389_ID
      load Men_389_BoundOpt.mat;
    case Net_913_ID
      load Men_913_BoundOpt.mat;
  end
  if RunType~=BOUND_OPT_SS
    if NetTypeID==Net_546_Meas_ID || NetTypeID==Net_389_ID || NetTypeID==Net_913_ID
      Boundary(:,3)=FinalBoundFlow;
    end
  end
end
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);

SegName=DataArray(:,1);
From=DataArray(:,3);
To=DataArray(:,4);
Len=DataArray(:,6);
SegType=DataArray(:,2);
VesNum=length(SegName);

%%%% 血管壁厚和粘滞度
WallTh=zeros(VesNum,1);
if NetTypeID~=Net_546_ID && NetTypeID~=Net_122_ID
  % 不是自适应后的546网络或者122网络，
  % 则管径从DataArray中读取，壁厚根据血管类型设置
  Diam=DataArray(:,5);
  for i=1:VesNum
%     if SegType(i)==1
%       % 不同类型血管的壁厚管径比根据546网络自适应后的结果得到
%       WallTh(i)=Diam(i)*0.2662;
%     elseif SegType(i)==2
%       WallTh(i)=Diam(i)*0.1812;
%     elseif SegType(i)==3
%       WallTh(i)=Diam(i)*0.0686;
%     else
%       WallTh(i)=Diam(i)*0.2;
%     end
    WallTh(i)=4;  % 为自适应模型设置
  end
  % 初始粘滞度统一设置为2mPas
  Visc=2*ones(VesNum,1);
else
  % 自适应血管网络，管径，壁厚，粘滞度从FuncPara中获取
  Diam=FuncPara(:,3);
  WallTh=FuncPara(:,20);
  Visc=FuncPara(:,9);
end

% 初始DampFactorName
DampFactorName='Init';
%%%% 粘滞度调整
DiamRatio=1.5;  % 管径调节的比例
switch DampFactor
  case Visc01
    Visc=Visc/10;
    DampFactorName='Visc01';
  case Visc10
    Visc=Visc*10;
    DampFactorName='Visc10';
  case ViscA01
    Visc(SegType==1)=Visc(SegType==1)/10;
    DampFactorName='ViscA01';
  case ViscA10
    Visc(SegType==1)=Visc(SegType==1)*10;
    DampFactorName='ViscA10';
  case ViscC01
    Visc(SegType==2)=Visc(SegType==2)/10;
    DampFactorName='ViscC01';
  case ViscC10
    Visc(SegType==2)=Visc(SegType==2)*10;
    DampFactorName='ViscC10';
  case ViscV01
    Visc(SegType==3)=Visc(SegType==3)/10;
    DampFactorName='ViscV01';
  case ViscV10
    Visc(SegType==3)=Visc(SegType==3)*10;
    DampFactorName='ViscV10';
  case ViscBD01
    BoundFlow(BoundType==0)=BoundFlow(BoundType==0)/10;
    DampFactorName='ViscBD01';
  case ViscBD10
    BoundFlow(BoundType==0)=BoundFlow(BoundType==0)*10;
    DampFactorName='ViscBD10';
  case DiamN
    Diam=Diam*DiamRatio;
    DampFactorName='DiamN';
  case DiamAN
    Diam(SegType==1)=Diam(SegType==1)*DiamRatio;
    DampFactorName='DiamAN';
  case DiamCN
    Diam(SegType==2)=Diam(SegType==2)*DiamRatio;
    DampFactorName='DiamCN';
  case DiamVN
    Diam(SegType==3)=Diam(SegType==3)*DiamRatio;
    DampFactorName='DiamVN';
  case Visc5
    Visc=Visc*5;
    DampFactorName='Visc5';
  case Visc02
    Visc=Visc/5;
    DampFactorName='Visc02';
  case DampVisc
    Visc=Visc*ViscRatio;
    DampFactorName=['Visc' num2str(ViscRatio)];
end
Boundary(:,3)=BoundFlow;
% load OptBound.mat;
% Boundary(:,3)=OptBound;

% AdjMode 管径调节模式
% 0: 保持不变
% 1：线性放大，2：指数放大，3：对数放大
% 4：线性缩小，5：指数缩小，6：对数缩小
[Diam DiamRatio]=AdjustDiam(Diam,SegType,0);

%%%%% 杨氏模量
% 所有网络的原始数据都没有杨氏模量，需根据血管类型进行设置
E=zeros(VesNum,1);
for i=1:VesNum
  if SegType(i)==1
    if ModelParam(18)==0
      E(i)=3.5e5;
    else
      % 杨氏模量设置根据(Salotto et al. 1986)
      if Diam(i)>25
        E(i)=6.3e5;
      elseif Diam(i)>15 && Diam(i)<25
        E(i)=2.6e5;
      else
        E(i)=1.6e5;
      end
    end
  elseif SegType(i)==3
    E(i)=3.88e5;
  else
    E(i)=3.7e5;
  end
end
%%%% 设置不同DampFactor下的杨氏模量
switch DampFactor
  case E01
    E=E/10;
    DampFactorName='E01';
  case E10
    E=E*10;
    DampFactorName='E10';
  case E02
    E=E/5;
    DampFactorName='E02';
  case E5
    E=E*5;
    DampFactorName='E5';
  case EA01
    E(SegType==1)=E(SegType==1)/10;
    DampFactorName='EA01';
  case EA10
    E(SegType==1)=E(SegType==1)*10;
    DampFactorName='EA10';
  case EC01
    E(SegType==2)=E(SegType==2)/10;
    DampFactorName='EC01';
  case EC10
    E(SegType==2)=E(SegType==2)*10;
    DampFactorName='EC10';
  case EV01
    E(SegType==3)=E(SegType==3)/10;
    DampFactorName='EV01';
  case EV10
    E(SegType==3)=E(SegType==3)*10;
    DampFactorName='EV10';
  case DampE
    E=E*ERatio;
    DampFactorName=['E' num2str(ERatio)];
end
DatMatrix=[SegName From To Len Diam WallTh SegType Visc E];

%%%% 只有自适应后的网络有FuncPara
if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  save('FuncPara.mat','FuncPara');
end

TaperRate=0;  % 锥度
%%%% 脉动性级别 %%%%
% 如果输入波形的类型是9，则需手动选择脉动性级别（在下面的程序中），其他情况下脉动性级别设为1
if ModelParam(INPUT_WAVE)~=9
  PulsatLevel=1;
end