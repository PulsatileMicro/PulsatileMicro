%% ����1Dģ�������ļ�(.in)��������
clear;clc;close all;
%% ModelParam: ģ�ͷ���������� %%%%
% 1 - Measured/Adapted flag 0: use Adapted diameter, 1: use Measured diameter
% (��������Ϊֻ��546������������⡣�ֽ�546��������������������������д���)
% 2 - Non-dimensionalization flag 1: Nondimensionalization
% 3 - Viscosity update strategy 0: No viscosity update 1: vitro 2: ESL vivo 3: No ESL vivo
% 4 - Riemann solver flag 1: Use riemann solver for the BioFlux
% 5 - Binary output flag 0: Output text format 1: Output binary format
% 6 - Viscoelastic wall 0: No viscous wall 1: New viscous model 2: Old viscous model
% 7 - Show CFL flag 0: Hide CFL number 1: Show CFL number
% 8 - Show Step Lapse flag 0: Hide Step Lapse 1: Show Step Lapse
% 9 - ODE solver flag 0: 1D_Exp 1: 1D_Imp 2: RLC_Exp 3: RLC_Imp 4: Steady
% State 5,6: Womersley (Not Ready) 7: RC_Exp 8: RC_Imp
% 10 - Input Waveform 0: Sublingual 1: Cat mesentery 2: Half sine 3:Impulse 4: Increasing Ramp
%   5: Decreasing Ramp 6: Sinusoidal 7: Multi-freq sine 8: Square 9: Sublingual Vel on Route A
%   10: Hypertension 0D-1D Coupling Input 11: Egg818 measured velocity
% 11 - Number of periods
% 12 - dt
% 13 - Period
% 14 - ODE Solver Rel tol
% 15 - ODE Solver Abs tol
% 16 - BoundType 0:u(mm/s) 1:q(ml/s) 2:p(mmHg)
% 17 - Introduce random phase difference to secondary boundaries (in 546 network)
% 18 - Linearity of Pressure-Area relationship: 0: Linear 1: Non-Linear
% 19 - Velocity profile: 0: Parabolic 1: Measured
% 20 - BG_Len_ratio
% 21 - BG_Mass_ratio
ModelParam=[0 0 0 1 0 0 0 1 4 0 4 0 0 1e-6 1e-6 1 0 0 1 10 100 1 0]; % 0DѪ������������
% ModelParam=[1 0 0 1 0 0 0 1 1 0 6 0 0 1e-6 1e-6 1 0 0 1 1e-2 1e-7]; % 1DѪ������������
% ModelParam=[1 0 0 1 0 0 0 1 0 0 4 0 0 1e-6 1e-6 0 0 0 1 10 100]; % �����ṹ�������

%%%% ���к궨�� %%%%
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID SymNet_ID Junc_ID Single_ID Tree_ID Salotto_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01
Macro
global ADAP NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE N_PERIOD DT PERIOD ODEREL ODEABS IN_BOUND_TYPE RANDOM_PHASE LIN_PA VEL_PROFILE BG_LEN_RATIO BG_MASS_RATIO

%%%% ���õ�ǰ������������ %%%%
NetTypeID=SymNet_ID;
if NetTypeID<100
  %%%% Damping�����ų�ģʽ: DampFactor
  DampFactor=DampInit;
end
% NetTypeName���������ļ�
[NetTypeName, ModelParam]=GetNetTypeName(NetTypeID,ModelParam);

%%%% ���沽�� %%%%
[dt, ModelParam]=GetDt(ModelParam);

%%%% 0D�����е������� %%%%
len_ratio=ModelParam(20);
mass_ratio=ModelParam(21);

%%%% �����ź����� %%%%
%%%% ���ݲ�ͬѪ�����͵����ź����� %%%%
ModelParam=GetPeriod(NetTypeID,ModelParam,1);
Period=ModelParam(PERIOD);

%%%% �����Լ��� %%%%
% �������������9�������ֶ�ѡ�������Լ���������ĳ����У������������Լ�����Ϊ1
if ModelParam(10)~=9
  PulsatLevel=1;
end

%% �����������Ѫ�����磬��ȡ/���ö�Ӧ��Ѫ������ %%%%
ModelParam(3)=0;    % ������Viscosity
switch NetTypeID
  case SymNet_ID
    InputDiam=100;
    order=8;    % orderָ����㵽ëϸѪ�ܲ�Ľ������ܽ���Ϊ2*order-1
    curOrder=0; cnt=0; lastNum=0;
    VesNum=2^order-1 + 2^(order-1)-1;
    From=zeros(VesNum,1);To=From;
    OrderInd=zeros(2*order-1,2);
    StartSeg=1;
    for i=1:2*order-1
      if i<=order
        OrderInd(i,1)=StartSeg;
        OrderInd(i,2)=OrderInd(i,1)+2^(i-1)-1;
        StartSeg=OrderInd(i,2)+1;
      else
        OrderInd(i,1)=StartSeg;
        OrderInd(i,2)=OrderInd(i,1)+2^(2*order-1-i)-1;
        StartSeg=OrderInd(i,2)+1;
      end
    end
    
%     Visc=2*ones(VesNum,1);
    VesType=4*ones(VesNum,1);
%     VesType(round(VesNum/2):end)=3;
%     VesType(OrderInd(6,1):OrderInd(6,2))=2;
    
    ratio_L=25;
    Visc=3;             %mPas
    ratio_D_L=0.75;
    ratio_D_R=0.75;
    ratio_D_Conv=0.75;
  case Tree_ID
    InputDiam=30;
    order=4;    % orderָ����㵽ëϸѪ�ܲ�Ľ���
    curOrder=0; cnt=0; lastNum=0;
    VesNum=2^order-1;
    OrderInd=zeros(2*order-1,2);
    StartSeg=1;
    for i=1:order
      OrderInd(i,1)=StartSeg;
      OrderInd(i,2)=OrderInd(i,1)+2^(i-1)-1;
      StartSeg=OrderInd(i,2)+1;
    end
    
    Visc=3*ones(VesNum,1);
    VesType=4*ones(VesNum,1);
%     VesType(OrderInd(order,1):OrderInd(order,2))=2;
    
%     ratio_L=25;
%     ratio_D_L=0.9182;
%     ratio_D_R=0.5876;
    ratio_L=25;
    ratio_D_L=0.75;
    ratio_D_R=0.75;

    %% �ӺϽṹ %
  case Junc_ID
    %       order=6;
    %       VesNum=2*order-1;
    VesNum=2;
    DiamRatio=1;
    InputDiam=40;               % um
    Diam=InputDiam*ones(VesNum,1);
    %       ratio_L_order=[30 40 50 60 70 80 70 60 50 40 30];
    ratio_L_order=50*ones(1,VesNum);
    VesType(1:order-1)=1;
    VesType(order)=2;
    VesType(order+1:VesNum)=3;
    VesType=4*ones(VesNum,1);
    for i=1:VesNum
      %         if i==1
      %           Diam(i)=InputDiam;
      %         elseif i<=order
      %           Diam(i)=Diam(i-1)*DiamRatio;
      %         else
      %           Diam(i)=Diam(i-1)/DiamRatio;
      %         end
      WallTh(i)=Diam(i)*0.2;
      %         if VesType(i)==1
      %           WallTh(i)=Diam(i)*0.2662;
      % %           WallTh(i)=WallTh(i);
      %         elseif VesType(i)==2
      %           WallTh(i)=Diam(i)*0.1812;
      % %           WallTh(i)=WallTh(i);
      %         elseif VesType(i)==3
      %           WallTh(i)=Diam(i)*0.0686;
      % %           WallTh(i)=WallTh(i);
      %         else
      % %           WallTh(i)=Diam(i)*0.2;
      %         end
    end
    Len=ratio_L_order'.*Diam;
    SegName=1:VesNum;
    %       Visc=[80 40 20 10 5 2.5 2.5 2.5 2.5 2.5 2];
    %       Visc(4:end)=Visc(4:end)/5;
    %       Visc(2:end)=Visc(2:end)/2;
    
    Visc=5*ones(VesNum,1);
    %       Visc(1)=Visc(1)*5;
    %       Visc(2)=Visc(2)*5;
    Hd=0.45;
  case Single_ID
    Len=10e3;                 % Length(um)       ճ����ЧӦ�볤���й�
    InputDiam=40;             % Diameter(um)
    Diam=InputDiam;
    WallTh=0.1*InputDiam;         % Wall thickness(um)
    Visc=2;                   % mPas
    VesType=1;
    
    VesNum=1;
    SegName=1;
    Hd=0.45;
end
% �����ṹ�ķ����������
RefCoeff=0.0;
ERatio=1;
ViscRatio=1;
TaperRate=0;
NumHisPt=3;
Freq=1;

switch ModelParam(16)
  case 0
    InputData=1;    % mm/s
    Vel=InputData;
  case 1
    Vel=1;
    InputData=Vel*pi/4*(InputDiam/1e3)^2; % mm^3/s
  case 2
    InputData=10;    % mmHg
    Vel=60;
end

%% %% �߽�����洢�����ʼ�� %%%%
% �˴��߽�ָ��������Ѫ�ܶεı߽磬���ǽ���ָѪ������ĳ���߽�
% ����Ѫ�ܶ��������ڴ��ı߽����ͣ���СΪVesNum*2
BCTypeAll=[];
% ���ڷֲ桢��ۡ��Ӻϱ߽磬BifurΪ����Ѫ�ܱ��
% ��������߽磬BifurΪ���������ͱ�ţ���u 0, u 1��0��1
% ��������߽磬Bifur��������
Bifur=zeros(VesNum,4);
% ���ڷֲ桢��ۡ��Ӻϱ߽磬BCVal��������
% ��������߽磬BCValΪ����߽��ֵ����a=PI*1e-4
% ��������߽磬BCValΪ����߽�R��T��ֵ
BCVal=zeros(VesNum,1);
% �߽����پ���
BoundData=[];

%% %% Ϊ����Ѫ�����ñ߽� %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% �Գ����� & Ѫ���� %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if NetTypeID==SymNet_ID || NetTypeID==Tree_ID
  for i=1:VesNum
    BCType=[];
    VesNumVec(i)=i;
    SegName(i)=i;
    if mod(i,2^(curOrder))==0
      curOrder=curOrder+1;
    end
    if curOrder<=order
      if i==1
        Diam(1)=InputDiam;
        Len(1)=InputDiam*ratio_L(1);
        Bifur(1,1:2)=[0 2];
        BoundData=[BoundData;InputData 1 0];
      elseif mod(i,2)==0
        Diam(i)=Diam(i/2)*ratio_D_L;
        Len(i)=Diam(i)*ratio_L;
        Bifur(i,1)=i/2;
        Bifur(i,2)=i+1;
        Bifur(i/2,3)=i;
        Bifur(i/2,4)=i+1;
      else
        Diam(i)=Diam((i-1)/2)*ratio_D_R;
        Len(i)=Diam(i)*ratio_L;
        Bifur(i,1)=(i-1)/2;
        Bifur(i,2)=i-1;
      end
      if curOrder==order
        if NetTypeID==SymNet_ID
          BCType=[BCType 'B' 'C'];
        else
%           BCType=[BCType 'B' 'R'];
%           BCVal(i)=30*133/(InputData*1e-9/256);
%           BCType=[BCType 'B' 'q'];
%           BCVal(i)=InputData*1e-9/(2^(order-1));
          BCType=[BCType 'B' 'R'];
          BCVal(i)=1e-5;
%           BCType=[BCType 'B' 'W'];
%           Bifur(i,3:4)=[30*133/(InputData*1e-9/256)/1e20 30*133/(InputData*1e-9/256)];
        end
      elseif i==1
        switch ModelParam(16)
          case 0
            BCType=[BCType 'u' 'B'];
            Bifur(i,1:2)=[0.45 2];
          case 1
            BCType=[BCType 'q' 'B'];
            Bifur(i,1:2)=[0.45 2];
            BCVal(i)=InputData*1e-9;
          case 2
            BCType=[BCType 'p' 'B'];
            Bifur(i,1:2)=[0.45 2];
        end
      else
        BCType=[BCType 'B' 'B'];
      end
    else
      if cnt==2^(2*order-curOrder-1)
        cnt=0;
        curOrder=curOrder+1;
      end
      lastNum=2*2^(2*order-curOrder-1);
      Diam(i)=Diam(i-lastNum+cnt)/ratio_D_Conv;
      Len(i)=Diam(i)*ratio_L;
      Bifur(i,1)=i-lastNum+cnt;
      Bifur(i,2)=i-lastNum+cnt+1;
      Bifur(i-lastNum+cnt,3)=i;
      Bifur(i-lastNum+cnt,4)=i-lastNum+cnt+1;
      Bifur(i-lastNum+cnt+1,3)=i-lastNum+cnt;
      Bifur(i-lastNum+cnt+1,4)=i;
      cnt=cnt+1;
      if i==VesNum
%         BCType=[BCType 'C' 'R'];
        BCType=[BCType 'C' 'W'];
%         BCVal(i)=RefCoeff;
%         BCVal(i)=30*133/(InputData*1e-9);
        Bifur(i,3:4)=[30*133/(InputData*1e-9)/1e20 30*133/(InputData*1e-9)];
      else
        BCType=[BCType 'C' 'C'];
      end
    end
    BCTypeAll=[BCTypeAll;BCType];
    Hd(i)=0.45;
  end
  for i=1:VesNum
    if VesType(i)==1
      WallTh(i)=Diam(i)*0.2662;
    elseif VesType(i)==2
      WallTh(i)=Diam(i)*0.1812;
    elseif VesType(i)==3
      WallTh(i)=Diam(i)*0.0686;
    else
      WallTh(i)=Diam(i)*0.1;
    end
  end
  Diam(VesType==1)=Diam(VesType==1);
  Diam(VesType==2)=Diam(VesType==2);
  Diam(VesType==3)=Diam(VesType==3);
  Len(VesType==1)=Len(VesType==1);
  Len(VesType==2)=Len(VesType==2);
  Len(VesType==3)=Len(VesType==3);
  %   Diam=Diam*2;
  %   Len=Len;
  
  %%%%%%%%%%%%%%%%%
  %%%% �Ӻ�Ѫ�� %%%%
  %%%%%%%%%%%%%%%%%
elseif NetTypeID==Junc_ID
  for i=1:VesNum
    if i==1
      %       Diam(i)=InputDiam;
      Bifur=[0 2 2 2];
      switch ModelParam(16)
        case 0
          BCTypeAll=['u' 'J'];
        case 1
          BCTypeAll=['q' 'J'];
        case 2
          BCTypeAll=['p' 'J'];
      end
    else
      %       Diam(i)=Diam(i-1)*DiamRatio(i);
      if i==VesNum
        Bifur=[Bifur;i-1 i-1 0 0];
        BCTypeAll=[BCTypeAll;'J' 'T'];
        BCVal(i)=RefCoeff;
      else
        Bifur=[Bifur;i-1 i-1 i+1 i+1];
        BCTypeAll=[BCTypeAll;'J' 'J'];
      end
    end
  end
  BoundData=[BoundData;InputData 1 0];
  %%%%%%%%%%%%%%%%%
  %%%% ����Ѫ�� %%%%
  %%%%%%%%%%%%%%%%%
elseif NetTypeID==Single_ID
  BCTypeAll=[];
  BCType=[];
  switch ModelParam(16)
    case 0
      BCType=['u' 'T'];
    case 1
      BCType=['q' 'T'];
    case 2
      BCType=['p' 'T'];
  end
  Bifur=[0 2];
  BCVal=RefCoeff;
  BCTypeAll=[BCTypeAll;BCType];
  BoundData=[BoundData;InputData 1 0];
end

%% %% ����Ѫ�ܲ��� %%%%
VesParam=zeros(22,VesNum);
VesParam(1,:)=Len'*1e-6;            % Length
VesParam(2,:)=Diam'*1e-6;           % Diameter
VesParam(3,:)=WallTh'*1e-6;         % Wall thickness
% ����Ѫ���������ò�ͬ������ģ��
for i=1:VesNum
  if NetTypeID>=100    % ����Ƿ���ʵѪ�����磬������ͳһ������ģ��
    if VesType(i)==1
      VesParam(4,i)=3.5e5;
    elseif VesType(i)==2
      VesParam(4,i)=3.5e5;
    else
%       VesParam(4,i)=3.5e5;
      VesParam(4,i)=5*Eh_Olufsen(Diam(i)/2/1e3/10)/10;
    end
  end
end
VesParam(4,:)=VesParam(4,:)*ERatio;
VesParam(5,:)=Eval_Phi(VesParam(4,:),Period);                  % Viscous wall modulus
if ModelParam(19)==0                % Alpha
  VesParam(6,:)=4/3;                % �ȼ���0Dģ�͵�����
elseif ModelParam(19)==1
  VesParam(6,:)=1.26;
end
VesParam(7,:)=Visc'*1e-3*ViscRatio; % Viscosity
if ModelParam(9)==2 || ModelParam(9)==3 || ModelParam(9)==7 || ModelParam(9)==8 % 0Dģ�ͣ���������
  VesParam(1:3,:)=VesParam(1:3,:)*1e3*len_ratio;
  VesParam([4 7],:)=VesParam([4 7],:)*mass_ratio;
end
VesParam(8,:)=Hd;                   % Discharge hematocrit
% ����Ǵ�Ѫ�����磬������SegName, From, To, Vel
if exist('From','var')
  VesParam(9,:)=SegName;            % Segment Name
  VesParam(10,:)=From;              % From nodes
  VesParam(11,:)=To;                % To nodes
  VesParam(12,:)=Vel;
end
% ����ǵ���Ѫ�ܣ����������٣�������Լ������ʱ����ʼ�׶���������
if NetTypeID==Single_ID
  VesParam(12,:)=Vel;
end
% ����1Dģ�͵�ճ����Ѫ�ܱ�ģ�Ͳ���
if ModelParam(6)==1
  VesParam(13,:)=Eval_Gamma(VesParam(2,:),VesParam(5,:),VesParam(3,:));  % Visc part of viscoelasticity
elseif ModelParam(6)==2
  VesParam(14,:)=Eval_GammaII(VesParam(2,:),VesParam(5,:));  % Visc part of viscoelasticity
else
  VesParam(13,:)=0;
  VesParam(14,:)=0;
end
if VesNum>10
  VesParam(15,:)=3;                   % q
  VesParam(16,:)=3;                   % L
else
  VesParam(15,:)=20;                   % q
  VesParam(16,:)=20;                   % L
end
% ȥ���ٻ�����
if ModelParam(2)
  VesParam(17,:)=1e-1;              % scale_lamda
  VesParam(18,:)=1e2;               % scale_u0
  VesParam(19,:)=1e-2;              % scale_r0
else
  VesParam(17,:)=1;
  VesParam(18,:)=1;
  VesParam(19,:)=1;
end
VesParam(20,:)=NumHisPt;            % Number of history points
VesParam(21,:)=TaperRate;           % Taper Rate

%% %% ���ɱ߽����������ļ� %%%%
% ģ�ʹ������ɵ��ļ��ж�ȡ�߽粨��
% Steady Stateģ�Ͳ���Ҫ
if ModelParam(9)~=4
  inFileName=[NetTypeName '_IN'];
  outFileName=[NetTypeName '_IN'];
  [i_artinput i_veninput]=GenBoundInput(inFileName,outFileName,BoundData,ModelParam,VesParam,Freq,PulsatLevel);
end

%% %% ����1Dģ�ͷ��������ȡ��.in�ļ� %%%%
save('ModelParam.mat','ModelParam');
save('VesParam.mat','VesParam');
save('NetType.mat','NetTypeName','NetTypeID','SymNet_ID','Single_ID','Tree_ID','Junc_ID');
save('Structure.mat','Diam');
save('Order.mat','order');
fileName=GenInputFile(NetTypeName,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);
system(['copy ' NetTypeName '*.* E:\Projects\HemoSim\HemoSim\']);
system(['copy ' NetTypeName '*.* E:\1DWin\' NetTypeName '\']);
system(['copy ModelParam.mat E:\1DWin\' NetTypeName '\']);
system(['copy VesParam.mat E:\1DWin\' NetTypeName '\']);
system(['copy NetType.mat E:\1DWin\' NetTypeName '\']);
system(['copy Structure.mat E:\1DWin\' NetTypeName '\']);
if NetTypeID<100
  % ����Ŀ¼���й�����ļ�
  system('del *.bcs *.in *.his');
end
if NetTypeID==Single_ID && ModelParam(10)==6
  % ���浥��Ѫ��ʱ����
  %   mkdir(['E:\1DWin\Single_D' num2str(Diam) '_f' num2str(Freq)]);
  %   system(['copy ' NetTypeName '*.* E:\1DWin\Single_D' num2str(Diam) '_f' num2str(Freq)]);
  %   system(['copy run_single.bat ' 'E:\1DWin\Single_D' num2str(Diam) '_f' num2str(Freq)]);
  %   mkdir(['E:\1DWin\Single_Visc' num2str(ViscRatio) '_f' num2str(Freq)]);
  %   system(['copy ' NetTypeName '*.* E:\1DWin\Single_Visc' num2str(ViscRatio) '_f' num2str(Freq)]);
  %   system(['copy run_single.bat ' 'E:\1DWin\Single_Visc' num2str(ViscRatio) '_f' num2str(Freq)]);
  mkdir(['E:\1DWin\Single_E' num2str(ERatio) '_f' num2str(Freq)]);
  system(['copy ' NetTypeName '*.* E:\1DWin\Single_E' num2str(ERatio) '_f' num2str(Freq)]);
  system(['copy run_single_HemoSim.bat ' 'E:\1DWin\Single_E' num2str(ERatio) '_f' num2str(Freq)]);
elseif NetTypeID==SymNet_ID && ModelParam(10)==6 % �����ź�����
  mkdir(['E:\1DWin\SymNet_f' num2str(Freq)]);
  system(['copy ' NetTypeName '*.* E:\1DWin\SymNet_f' num2str(Freq)]);
  system(['copy run_sym_HemoSim.bat ' 'E:\1DWin\SymNet_f' num2str(Freq)]);
end


