clear;clc;close all;warning off;
%% ΢ѭ��ģ�ͳ���������
%%%% ���к궨�� %%%%
Macro   % ����궨����ļ�
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global SIM_PREP FREQ_ANAL BOUND_OPT STRUCT_ADAP
%% �ó��������������ģʽ
% SIM_PREP - C���Ե�0D,1D,��̬ģ�ͷ����׼��(��������)
% FREQ_ANAL - Ƶ�����
% BOUND_OPT - �߽��Ż�
% STRUCT_ADAP - ����Ӧģ��
RunType=STRUCT_ADAP;
%%%% ����ѡ�� %%%%
NetTypeID=Net_122_ID;
%%%% ����NetTypeID��ȡ��Ӧ��Ѫ����������
[DataArray,Boundary,FuncPara,DatMatrix]=LoadRealNetData(RunType,NetTypeID);
%%%% ������������������Ӧ����
switch RunType
  case SIM_PREP
    [SS_Press,SS_Flow]=SIM_PREP_FUNC(NetTypeID,Boundary,DatMatrix);
  case FREQ_ANAL
    FREQ_ANAL_FUNC();
  case BOUND_OPT
    [OptBound,InvSeg]=BOUND_OPT_FUNC(NetTypeID,Boundary,DatMatrix);
  case STRUCT_ADAP
    [AdapCoeff,FVal,OptDiam,OptWallTh,Sm,Sc]=...
      STRUCT_ADAP_FUNC(NetTypeID,DataArray,Boundary,FuncPara,DatMatrix);
end

% ���ѳ������н���
load chirp
sound(y,Fs)