%% ���ݷ����Ѫ���������ͣ���������
%%%% �����Ҫ���ĵ�������ݣ����޸�������� %%%%
function [DataArray,Boundary,FuncPara,DatMatrix]=LoadRealNetData(RunType,NetTypeID)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global SIM_PREP FREQ_ANAL BOUND_OPT STRUCT_ADAP

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
[DataArray,Boundary,FuncPara]=ReadData(DatFile,PrnFile);
SegName=DataArray(:,1);
From=DataArray(:,3);
To=DataArray(:,4);
Len=DataArray(:,6);
SegType=DataArray(:,2);
VesNum=length(SegType);
% ʹ���Ż��õ���Ѫ������߽�
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

%%%% Ѫ�ܱں��ճ�Ͷ�
WallTh=zeros(VesNum,1);
if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  % ������Ӧ���Ѫ�����磬�ܾ����ں�ճ�Ͷȴ�FuncPara�л�ȡ
  Diam=FuncPara(:,3);     % ����Ӧ��Ѫ��
  WallTh=FuncPara(:,20);
  Visc=FuncPara(:,9);
elseif NetTypeID==Egg_818_ID || NetTypeID==Egg_636_ID
  Diam=DataArray(:,5);
  WallTh=Diam/2*0.02;     % ����Pries 2005����
  % ��ʼճ�Ͷ�ͳһ����Ϊ2mPas
  Visc=2*ones(VesNum,1);
else
  % ��������Ӧ���546�������122���磬
  % ��ܾ���DataArray�ж�ȡ���ں����Ѫ����������
  Diam=DataArray(:,5);
  for i=1:VesNum
    if RunType==STRUCT_ADAP
      WallTh(i)=4;  % Ϊ����Ӧģ������ͳһ�ĳ�ֵ����λum
    else
      if SegType(i)==1
        % ��ͬ����Ѫ�ܵıں�ܾ��ȸ���546��������Ӧ��Ľ���õ�
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
  % ��ʼճ�Ͷ�ͳһ����Ϊ2mPas
  Visc=2*ones(VesNum,1);
end


%%%%% ����ģ��
% ���������ԭʼ���ݶ�û������ģ���������Ѫ�����ͽ�������
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
        %       % ����ģ�����ø���(Salotto et al. 1986)
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