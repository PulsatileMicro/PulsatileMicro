clear;clc;close all;warning off;dbstop if error;
%% 计时开始
SimulateInitT = clock;
%% 微循环模型程序主函数
%%%% 运行宏定义 %%%%
Macro   % 定义宏定义的文件
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global SIM_PREP FREQ_ANAL BOUND_OPT STRUCT_ADAP
%% 该程序具有以下运行模式
% SIM_PREP - C语言的0D,1D,稳态模型仿真的准备(数据生成)
% FREQ_ANAL - 频域分析
% BOUND_OPT - 边界优化
% STRUCT_ADAP - 自适应模型
RunType=STRUCT_ADAP;
%%%% 网络选择 %%%%
NetTypeID=Net_122_ID;
%%%% 根据NetTypeID读取相应的血管网络数据
[DataArray,Boundary,FuncPara,DatMatrix]=LoadRealNetData(RunType,NetTypeID);
%%%% 根据运行类型运行相应函数
switch RunType
  case SIM_PREP
    [SS_Press,SS_Flow]=SIM_PREP_FUNC(NetTypeID,Boundary,DatMatrix);
  case FREQ_ANAL
    FREQ_ANAL_FUNC();
  case BOUND_OPT
    [OptBound,InvSeg]=BOUND_OPT_FUNC(NetTypeID,Boundary,DatMatrix);
  case STRUCT_ADAP
    [AdapCoeff,FVal,OptDiam,OptWallTh,PO2,MeanP,tau,O,Visc,Sm,Sc]=...
      STRUCT_ADAP_FUNC(NetTypeID,DataArray,Boundary,FuncPara,DatMatrix);
end

%% 计时结束
disp(['toc：总共运行时间：', num2str(etime(clock, SimulateInitT))]);
% 提醒程序运行结束
load chirp
sound(y,Fs)