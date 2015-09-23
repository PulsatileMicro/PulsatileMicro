%% ���ݷ����Ѫ���������ͣ���������
%%%% �����Ҫ���ĵ�������ݣ����޸�������� %%%%
function [ModelParam,DataArray,Boundary,FuncPara,DatMatrix,VesNum,TaperRate,PulsatLevel,DampFactorName]...
=LoadRealNetData(NetTypeID,ModelParam,DampFactor,RunType,ViscRatio,ERatio)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN E5 E02 Visc5 Visc02 DampE DampVisc
global ADAP NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE N_PERIOD DT PERIOD ODEREL ODEABS IN_BOUND_TYPE RANDOM_PHASE LIN_PA VEL_PROFILE BG_LEN_RATIO BG_MASS_RATIO USE_OPTBOUND
global ONED_PREP FREQ_ANAL BOUND_OPT_SS BOUND_OPT_PULSATILE

switch NetTypeID
    %% 546����Ӧ������
  case Net_546_ID
    DatFile='Men_546.dat';
    PrnFile='Men_546.prn';
    %% 546����Ӧǰ����
  case Net_546_Meas_ID
    DatFile='Men_546.dat';
    PrnFile='';
    %% 389Ѫ������
  case Net_389_ID
    DatFile='Men_389.dat';    % Network data file name
    PrnFile='';
    %% 913Ѫ������
  case Net_913_ID
    DatFile='Men_913.dat';    % Network data file name
    PrnFile='';
    %% ����818Ѫ������
  case Egg_818_ID
    DatFile='Vit_818.dat';
    PrnFile='';
    %% ����636Ѫ������
  case Egg_636_ID
    DatFile='Vit_636.dat';
    PrnFile='';
    %% CAMѪ������
  case Egg_CAM_ID
    DatFile='CAM_7128.DAT';    % Network data file name
    PrnFile='';
    %% SubCAMѪ������
  case Sub_CAM_ID
    DatFile='SubCAM.DAT';    % Network data file name
    PrnFile='';
    %% 122����Ѫ������
  case Net_122_ID
    DatFile='Men_122.dat';
    PrnFile='Men_122.prn';
end

%% ����������ʵѪ�����磬���ļ���ȡ����
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
% ����CAM���磬������߽�����ֵ
% ��Ϊ�߽�����ֵ����ʵ��ģ�������������-�ܾ��ľ����ϵ�����ģ���Ҫ����΢��
% if NetTypeID==Egg_CAM_ID
%   Boundary(10:11,3)=0.35*Boundary(10:11,3);
% end

% �߽�����
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

%%%% Ѫ�ܱں��ճ�Ͷ�
WallTh=zeros(VesNum,1);
if NetTypeID~=Net_546_ID && NetTypeID~=Net_122_ID
  % ��������Ӧ���546�������122���磬
  % ��ܾ���DataArray�ж�ȡ���ں����Ѫ����������
  Diam=DataArray(:,5);
  for i=1:VesNum
%     if SegType(i)==1
%       % ��ͬ����Ѫ�ܵıں�ܾ��ȸ���546��������Ӧ��Ľ���õ�
%       WallTh(i)=Diam(i)*0.2662;
%     elseif SegType(i)==2
%       WallTh(i)=Diam(i)*0.1812;
%     elseif SegType(i)==3
%       WallTh(i)=Diam(i)*0.0686;
%     else
%       WallTh(i)=Diam(i)*0.2;
%     end
    WallTh(i)=4;  % Ϊ����Ӧģ������
  end
  % ��ʼճ�Ͷ�ͳһ����Ϊ2mPas
  Visc=2*ones(VesNum,1);
else
  % ����ӦѪ�����磬�ܾ����ں�ճ�Ͷȴ�FuncPara�л�ȡ
  Diam=FuncPara(:,3);
  WallTh=FuncPara(:,20);
  Visc=FuncPara(:,9);
end

% ��ʼDampFactorName
DampFactorName='Init';
%%%% ճ�Ͷȵ���
DiamRatio=1.5;  % �ܾ����ڵı���
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

% AdjMode �ܾ�����ģʽ
% 0: ���ֲ���
% 1�����ԷŴ�2��ָ���Ŵ�3�������Ŵ�
% 4��������С��5��ָ����С��6��������С
[Diam DiamRatio]=AdjustDiam(Diam,SegType,0);

%%%%% ����ģ��
% ���������ԭʼ���ݶ�û������ģ���������Ѫ�����ͽ�������
E=zeros(VesNum,1);
for i=1:VesNum
  if SegType(i)==1
    if ModelParam(18)==0
      E(i)=3.5e5;
    else
      % ����ģ�����ø���(Salotto et al. 1986)
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
%%%% ���ò�ͬDampFactor�µ�����ģ��
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

%%%% ֻ������Ӧ���������FuncPara
if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  save('FuncPara.mat','FuncPara');
end

TaperRate=0;  % ׶��
%%%% �����Լ��� %%%%
% ������벨�ε�������9�������ֶ�ѡ�������Լ���������ĳ����У�����������������Լ�����Ϊ1
if ModelParam(INPUT_WAVE)~=9
  PulsatLevel=1;
end