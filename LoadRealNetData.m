%% 根据仿真的血管网络类型，导入数据
%%%% 如果需要更改导入的数据，需修改这个函数 %%%%
function [DataArray,Boundary,FuncPara,DatMatrix]=LoadRealNetData(RunType,NetTypeID)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global SIM_PREP FREQ_ANAL BOUND_OPT STRUCT_ADAP

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
[DataArray,Boundary,FuncPara]=ReadData(DatFile,PrnFile);
SegName=DataArray(:,1);
From=DataArray(:,3);
To=DataArray(:,4);
Len=DataArray(:,6);
SegType=DataArray(:,2);
VesNum=length(SegType);
% 使用优化得到的血管网络边界
switch NetTypeID
  case Net_546_Meas_ID
    load Men_546_Meas_BoundOpt.mat;
  case Net_389_ID
    load Men_389_BoundOpt.mat;
  case Net_913_ID
    load Men_913_BoundOpt.mat;
  case Egg_818_ID
    load Egg_818_BoundOpt.mat;
  case Egg_636_ID
    load Egg_636_BoundOpt.mat;
end
if exist('FinalBoundFlow')
  Boundary(:,3)=FinalBoundFlow;
end

%%%% 血管壁厚和粘滞度
WallTh=zeros(VesNum,1);
if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  % 是自适应后的血管网络，管径，壁厚，粘滞度从FuncPara中获取
  Diam=FuncPara(:,3);     % 自适应后血管
  WallTh=FuncPara(:,20);
  Visc=FuncPara(:,9);
elseif NetTypeID==Egg_818_ID || NetTypeID==Egg_636_ID
  Diam=DataArray(:,5);
  WallTh=Diam/2*0.02;     % 根据Pries 2005估算
  % 初始粘滞度统一设置为2mPas
  Visc=2*ones(VesNum,1);
else
  % 不是自适应后的546网络或者122网络，
  % 则管径从DataArray中读取，壁厚根据血管类型设置
  Diam=DataArray(:,5);
  for i=1:VesNum
    if RunType==STRUCT_ADAP
      WallTh(i)=4;  % 为自适应模型设置统一的初值，单位um
    else
      if SegType(i)==1
        % 不同类型血管的壁厚管径比根据546网络自适应后的结果得到
        WallTh(i)=Diam(i)*0.2662;
      elseif SegType(i)==2
        WallTh(i)=Diam(i)*0.1812;
      elseif SegType(i)==3
        WallTh(i)=Diam(i)*0.0686;
      else
        WallTh(i)=Diam(i)*0.2;
      end
    end
  end
  % 初始粘滞度统一设置为2mPas
  Visc=2*ones(VesNum,1);
end


%%%%% 杨氏模量
% 所有网络的原始数据都没有杨氏模量，需根据血管类型进行设置
E=zeros(VesNum,1);
switch NetTypeID
  case {Egg_818_ID,Egg_636_ID,Egg_CAM_ID}
    for i=1:VesNum
      E(i)=Calc_E_From_Dist(3e-4,WallTh(i)/1e6,Diam(i)/1e6);
    end
  otherwise
    for i=1:VesNum
      if SegType(i)==1
        %     if ModelParam(18)==0
        %       E(i)=3.5e5;
        %     else
        %       % 杨氏模量设置根据(Salotto et al. 1986)
        %       if Diam(i)>25
        %         E(i)=6.3e5;
        %       elseif Diam(i)>15 && Diam(i)<25
        %         E(i)=2.6e5;
        %       else
        %         E(i)=1.6e5;
        %       end
        %     end
        E(i)=3.5e5;
      elseif SegType(i)==3
        E(i)=3.88e5;
      else
        E(i)=3.7e5;
      end
    end
end

DatMatrix=[SegName From To Len Diam WallTh SegType Visc E];