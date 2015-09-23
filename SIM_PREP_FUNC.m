function [SS_Press,SS_Flow] = SIM_PREP_FUNC(NetTypeID,Boundary,DatMatrix)
global DampInit DampVisc DampE
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP Sparse_SS
global NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE IN_BOUND_TYPE RANDOM_PHASE VEL_PROFILE USE_OPTBOUND BIFURLOSS
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global HUMAN CATMEN HALFSIN IMPULSE INCRAMP DECRAMP SINE MULTI_FREQ_SIN SQUARE HUMAN_ROUTEA_546 HYPER VIT

% �����������Macro.m
ModelParam(NONDIM)=0;
ModelParam(VISC_UPD)=0;
ModelParam(RIEM)=1;
ModelParam(BINOUT)=0;
ModelParam(VISCOELAS)=1;
ModelParam(CFL)=0;
ModelParam(STEP_LAPSE)=1;
ModelParam(ODESOLVER)=RC_IMP;
switch NetTypeID
  case {Egg_818_ID,Egg_636_ID,Egg_CAM_ID,Sub_CAM_ID}
    ModelParam(INPUT_WAVE)=VIT;
  otherwise
    ModelParam(INPUT_WAVE)=HUMAN;
end
ModelParam(IN_BOUND_TYPE)=1; % TODO:use macro name
ModelParam(RANDOM_PHASE)=0;
ModelParam(VEL_PROFILE)=1;
ModelParam(USE_OPTBOUND)=1;
ModelParam(BIFURLOSS)=1;

% ModelParam=[0 0 1 0 0 0 1 ONED_IMP 0 4 0 0 1e-6 1e-6 1 0 0 1 10 100 1 0]; % 0DѪ������������
%% Step 1: ׼������
%%%% 1.1 �����������ͣ��������������� %%%%
% TODO: ����ȡ��
% ����nCycle���Ķ�����
switch NetTypeID
  case {Egg_818_ID,Egg_636_ID,Egg_CAM_ID,Sub_CAM_ID}
    nCycle=12;    % Egg Period=280ms, in total 3360 steps
  otherwise
    nCycle=4;       % Mes Period=800ms, in total 3200 steps
end

%%%% 1.2 Damping�����ų�ģʽ: DampFactor %%%%
% DampFactor: DampInit, DampVisc, DampE
DampFactor=DampVisc;
ViscRatio=0.1;
ERatio=0.1;
DampPara=[DampFactor,ViscRatio,ERatio];
DampFactorName=GenDampFactorName(DampFactor,ViscRatio,ERatio);
fprintf('DampType=%s',DampFactorName);

%%%% 1.3 �����ź����ڱ��� %%%%
% ����=����/nFreqTimes����������
switch NetTypeID
  case {Egg_818_ID,Egg_636_ID,Egg_CAM_ID,Sub_CAM_ID}
    nFreqTimes=1;   % ����Egg��ȡ��ԭֵ
  otherwise
%     nFreqTimes=4;   % ����Mes����Ϊ�������벨�β������壬HRԼΪ75bpm��ȡ��4����300bpm��Ϊ����
    nFreqTimes=1;     % ����
end

% ���ݲ�ͬѪ�����ͣ���ͬ���ڱ����������ź����� %%%%
Period=GetPeriod(NetTypeID,ModelParam)/nFreqTimes;
nCycle=nCycle*nFreqTimes;

%%%% 1.4 �߽�ֵ������ %%%%
% ����߽�����Rʱ��ʵ���������ڱ߽�Ѫ�ܴ���õ��������������е���
switch NetTypeID
  case {Egg_818_ID,Egg_636_ID,Egg_CAM_ID,Sub_CAM_ID}
    OutBoundRatio=1;    % ����Egg����������Ҫ���ڣ���Ϊ��������ںܽӽ�������
  otherwise
%     OutBoundRatio=100;  % ����Mes����Ҫ���ڣ������ڶ����޷�ȷ��(�����ȷ���TODO)���ݶ�100
    OutBoundRatio=1;    % ����
end

%%%% 1.5 ���沽�� %%%%
dt=GetDt(ModelParam(ODESOLVER));

%%%% 1.6 0D�����е������� %%%%
if ModelParam(ODESOLVER)==RLC_EXP || ModelParam(ODESOLVER)==RLC_IMP ||...
    ModelParam(ODESOLVER)==RC_EXP || ModelParam(ODESOLVER)==RC_IMP
  len_ratio=10;
  mass_ratio=100;
else
  len_ratio=1;
  mass_ratio=1;
end

%%%% 1.7 NetTypeName��������.in�ļ� %%%%
% NetTypeName���ڸ�������ļ�����
NetTypeName=GetNetTypeName(NetTypeID);

%%%% 1.8 SolverName %%%%
SolverName=GetSolverName(ModelParam(ODESOLVER));

%% Step 2: ��̬Ѫ������
% ���ж�̬����ǰ���ȷ�����̬Ѫ������ѧ���ԣ�����������Ϊ��̬ģ�����ñ߽���������Ҫ������߽�����
[SS_Press,SS_Flow,SS_Vel,SS_DeltaP,SS_Visc,SS_Hd]=LinEqu(NetTypeID,DatMatrix,Boundary,DampPara);
VesNum=length(SS_Press);

%% Step 3: ����1Dģ�ͱ߽�����洢�����ʼ��
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
SegName=DatMatrix(:,1);
From=DatMatrix(:,2);
To=DatMatrix(:,3);
Len=DatMatrix(:,4);
Diam=DatMatrix(:,5);
WallTh=DatMatrix(:,6);
SegType=DatMatrix(:,7);
Visc=DatMatrix(:,8);
E=DatMatrix(:,9);

%%%% Ϊ����Ѫ�����ñ߽�
for i=1:VesNum
  BCType=[];  % ��ʱ�������洢�߽�����
  %%% ��������߽�
  inInd1=find(From(i)==From); % ��ѯ�����Ӹö�Ѫ����������Ѫ��
  inInd2=find(From(i)==To);   % ��ѯ����ö�Ѫ������Ѫ��
  if length(inInd1)==2
    % �����2��Ѫ�ܾ��и���㣬˵���ýڵ���һ���ֲ�ڵ�
    % ��ʱ�����ҽ���һ��Ѫ������õ㣬��inInd2
    BCType=[BCType 'B'];
    inInd1(inInd1==i)=[]; % ɾ����������Ѫ�ܵı��
    Bifur(i,1:2)=[inInd2 inInd1];
  elseif length(inInd2)==2
    % �����2��Ѫ���������㣬˵���ýڵ���һ����۽ڵ�
    BCType=[BCType 'C'];
    Bifur(i,1:2)=inInd2;
  elseif length(inInd2)==0
    % ���û��Ѫ���������㣬˵���ö�Ѫ����һ������߽�
    switch ModelParam(16)
      case 0
        BCType=[BCType 'u'];
        BCVal(i)=SS_Vel(i);  % Velocity in mm/s
        Bifur(i,1:2)=[SS_Hd(i) 2];
      case 1
        BCType=[BCType 'q'];
        BCVal(i)=SS_Vel(i);  % Flow in mm^3/s, �ܾ���GenBoundInput�������ٳ�
        if ModelParam(ODESOLVER)==RLC_EXP || ModelParam(ODESOLVER)==RLC_IMP ||...
            ModelParam(ODESOLVER)==RC_EXP || ModelParam(ODESOLVER)==RC_IMP   % �����0Dģ��
          BCVal(i)=BCVal(i);
        elseif ModelParam(ODESOLVER)==SS || ModelParam(ODESOLVER)==Sparse_SS
          % �������̬ģ�ͣ��߽�ֵ��Ϊ����ֵ
          BCVal(i)=SS_Flow(i)/60/1e12;
        end
        Bifur(i,1:2)=[SS_Hd(i) 2];
      case 2
        BCType=[BCType 'p'];
        BCVal(i)=SS_Press(i);
        Bifur(i,1:2)=[SS_Hd(i) 2];
    end
    
%     if NetTypeID==Egg_818_ID && (From(i)~=685 && From(i)~=682)
%       BoundData=[BoundData;BCVal(i) i 1];
%     elseif NetTypeID==Egg_CAM_ID && From(i)==1266  % Egg_CAM
%       BoundData=[BoundData;BCVal(i) i 1];
%       %       elseif NetTypeID==Egg_636_ID && VesType(i)==3
%     elseif NetTypeID==Egg_636_ID && From(i)~=137
%       BoundData=[BoundData;BCVal(i) i 1];
%       %       elseif NetTypeID==Net_546_ID && From(i)~=830
%       %         BoundData=[BoundData;BCVal(i) i 1];
%     else
%       BoundData=[BoundData;BCVal(i) i 0];
%     end
    BoundData=[BoundData;BCVal(i) i 0];
  else
    % �����������������Ϊ�Ӻ�Ѫ��
    BCType=[BCType 'J'];
    Bifur(i,1:2)=inInd2;
  end
  
  %%% ��������߽�
  outInd1=find(To(i)==From);  % ��ѯ�ö�Ѫ�������Ѫ��
  outInd2=find(To(i)==To);    % ��ѯ��������ö�Ѫ���յ��Ѫ��
  if length(outInd1)==2
    % ����ö�Ѫ����������Ѫ�ܣ�˵���ýڵ���һ���ֲ�ڵ�
    BCType=[BCType 'B'];
    Bifur(i,3:4)=outInd1;
  elseif length(outInd2)==2
    % �����2��Ѫ�ܵ�����յ㣬˵���ýڵ���һ����۽ڵ�
    BCType=[BCType 'C'];
    outInd2(outInd2==i)=[];
    Bifur(i,3:4)=[outInd2 outInd1];
  elseif length(outInd1)==0
    % ������յ㲻���κ�Ѫ�ܵ���㣬˵���ýڵ���һ�������߽�
    if ModelParam(ODESOLVER)~=SS && ModelParam(ODESOLVER)~=Sparse_SS % ������̬ģ�ͣ�����R�߽�
      BCType=[BCType 'R'];
      % ��ʵ��SS_FlowӦ�ø�������BCValӦ�ø�С
      % SS_Flow��������OutBoundRatio�仯
      BCVal(i)=(SS_Press(i)-SS_DeltaP(i)/2)/SS_Flow(i)/OutBoundRatio*133*60*1e12;
      if ModelParam(ODESOLVER)==RLC_EXP || ModelParam(ODESOLVER)==RLC_IMP ||...
          ModelParam(ODESOLVER)==RC_EXP || ModelParam(ODESOLVER)==RC_IMP
        BCVal(i)=BCVal(i)/1e9*mass_ratio/len_ratio^3;
      end
      % ����DampFactor���ڱ߽�����
      % ճ�Ͷȵ���ʱ���߽�����ֵҲ�������
      % ���ǲ��������Է��̷���ǰ���ڣ�����ᵼ����̬�붯̬ģ��Ѫѹˮƽ��һ��
      % �Ƿ�Ӧ�õ����������ۣ�����Ӧ�÷�Visc01��Visc01BD01����
      if DampFactor==DampVisc
        BCVal(i)=BCVal(i)*ViscRatio;
      end
    else
      % ��̬ģ��
      ind=find(To(i)==Boundary(:,1));
      if Boundary(ind,2)==0
        % ����߽�����Ϊѹ����������ѪѹΪ�߽�
        BCType=[BCType 'R'];
        BCVal(i)=SS_Press(i);
      else
        % ����߽�����Ϊ����������������Ϊ�߽磬����Ϊ��
        BCType=[BCType 'q'];
        BCVal(i)=-SS_Flow(i)/60/1e12;
      end
    end
  else
    % �����������������Ϊ�Ӻ�Ѫ��
    BCType=[BCType 'J'];
    Bifur(i,3:4)=outInd1;
  end
  BCTypeAll=[BCTypeAll;BCType];
end

%% 4. ����Ѫ�ܲ��� %%%%
VesParam=zeros(22,VesNum);
VesParam(1,:)=Len'*1e-6;            % Length
VesParam(2,:)=Diam'*1e-6;           % Diameter
VesParam(3,:)=WallTh'*1e-6;         % Wall thickness
if DampFactor==DampE
  VesParam(4,:)=E/ERatio;
else
  VesParam(4,:)=E;
end
VesParam(5,:)=Eval_Phi(VesParam(4,:),Period);% Viscous wall modulus
if ModelParam(19)==0                % Alpha
  VesParam(6,:)=4/3;                % alpha=4/3ʱ���ȼ���0Dģ�͵�����
elseif ModelParam(19)==1
  VesParam(6,:)=1.26;
end
VesParam(7,:)=SS_Visc'*1e-3; % Viscosity
if ModelParam(ODESOLVER)==RC_IMP || ModelParam(ODESOLVER)==RC_EXP ||...
    ModelParam(ODESOLVER)==RLC_IMP || ModelParam(ODESOLVER)==RLC_EXP % 0Dģ�ͣ���������
  VesParam(1:3,:)=VesParam(1:3,:)*1e3*len_ratio;
  VesParam([4 7],:)=VesParam([4 7],:)*mass_ratio;
end
VesParam(8,:)=SS_Hd;              % Discharge hematocrit
VesParam(9,:)=SegName;            % Segment Name
VesParam(10,:)=From;              % From nodes
VesParam(11,:)=To;                % To nodes
VesParam(12,:)=SS_Vel;
% ����1Dģ�͵�ճ����Ѫ�ܱ�ģ�Ͳ���
if ModelParam(VISCOELAS)==1
  VesParam(13,:)=Eval_Gamma_New(VesParam(2,:),VesParam(5,:));  % Visc part of viscoelasticity
elseif ModelParam(VISCOELAS)==2
  VesParam(14,:)=Eval_Gamma_Old(VesParam(2,:),VesParam(5,:),VesParam(3,:));  % Visc part of viscoelasticity
else
  VesParam(13,:)=0;
  VesParam(14,:)=0;
end
if VesNum>10  % ����Ѫ�����磬1Dģ���е�q,L������Ϊ3
  VesParam(15,:)=3;                   % q
  VesParam(16,:)=3;                   % L
else
  VesParam(15,:)=20;                   % q
  VesParam(16,:)=20;                   % L
end
% ȥ���ٻ�����
if ModelParam(NONDIM)
  VesParam(17,:)=1e-1;              % scale_lamda
  VesParam(18,:)=1e2;               % scale_u0
  VesParam(19,:)=1e-2;              % scale_r0
else
  VesParam(17,:)=1;
  VesParam(18,:)=1;
  VesParam(19,:)=1;
end
if ModelParam(ODESOLVER)==ONED_EXP || ModelParam(ODESOLVER)==ONED_IMP
  NumHisPt=3;
else
  NumHisPt=1;
end
VesParam(20,:)=NumHisPt;            % Number of history points
VesParam(22,:)=SegType;             % Ѫ�ܶ�����
VesParam(23,:)=ModelParam(VISCOELAS);

%% 5. ���ɱ߽����������ļ� %%%%
% ģ�ʹ������ɵ��ļ��ж�ȡ�߽粨��
% Steady Stateģ�Ͳ���Ҫ
if ModelParam(ODESOLVER)~=SS && ModelParam(ODESOLVER)~=Sparse_SS
  inFileName=[NetTypeName '_IN'];
  %   outFileName=[NetTypeName '_OUT'];
  outFileName=[NetTypeName '_IN'];
  [i_artinput i_veninput]=...
    GenBoundInput(inFileName,outFileName,BoundData,ModelParam,...
    VesParam,nFreqTimes,0,dt,Period,nCycle,len_ratio,mass_ratio);
end

%     % ȷ������Damping����£��߽�����һ�¡���Ҫ�ǳ��߽�Rһ��
%     switch NetTypeID
%       case Net_546_Meas_ID
%         load Men_546_Meas_BCVal.mat;
%       case Net_546_ID
%         load Men_546_BCVal.mat;
%       case Net_913_ID
%         load Men_913_BCVal.mat;
%       case Net_389_ID
%         load Men_389_BCVal.mat;
%     end
%     for i=1:VesNum
%       if BCVal(i)>500 % ����500��Ӧ���Ǳ߽�����ֵ��С��500��Ӧ���Ǳ߽�����ֵ
%         if ModelParam(ODESOLVER)==RLC_EXP || ModelParam(ODESOLVER)==RLC_IMP ||...
%             ModelParam(ODESOLVER)==RC_EXP || ModelParam(ODESOLVER)==RC_IMP
%           BCVal(i)=BCVal(i)/1e9*mass_ratio/len_ratio^3;
%         end
%       end
%     end

%% 6. ����1Dģ�ͷ��������ȡ��.in�ļ� %%%%
fileName=GenInputFile(NetTypeName,VesParam,BCTypeAll,Bifur,BCVal,ModelParam,DampFactorName,dt,nCycle,Period,len_ratio,mass_ratio);
% ������Ԥ�������Ϣ�������������ͬһĿ¼���Ա�Review�������������
save('ModelParam.mat','ModelParam','dt','Period','nCycle','len_ratio','mass_ratio');
save('VesParam.mat','VesParam');
save('NetType.mat','NetTypeName','NetTypeID','Net_546_ID','Net_546_Meas_ID','Egg_818_ID','Net_122_ID','Net_389_ID','Net_913_ID','Egg_CAM_ID','Sub_CAM_ID','Egg_636_ID');
save('Structure.mat','Diam');
save('SegName.mat','SegName');
% ��Ԥ������������ָ��λ�ý��з���
% ���ݷ����������Ŀ¼��������1.������, 2.���������, 3.Damp������, 4.�Ƿ�Viscoelastic wall 5.����, 6.�߽���ڱ���
FolderName=[NetTypeName '_' SolverName '_' DampFactorName '_Gamma' int2str(ModelParam(VISCOELAS)) '_bpm' int2str(60/Period)...
  '_BDFlow' int2str(OutBoundRatio) '_BifurLoss' int2str(ModelParam(BIFURLOSS))];
system(['mkdir E:\1DWin\' FolderName]);
system(['mkdir PulseAdapDIR']);
system(['copy ' NetTypeName '*.* E:\Projects\HemoSim\HemoSim\']);
system(['copy ' NetTypeName '*.* E:\Projects\HemoSim_Sundials2.6\HemoSim\']);
system(['copy ' NetTypeName '*.* E:\1DWin\' FolderName '\']);
system(['copy ' NetTypeName '*.* PulseAdapDIR\']);
% ����������Ŀ¼��
system(['copy LinEquDatFile.mat E:\1DWin\' FolderName '\']);
system(['copy ModelParam.mat E:\1DWin\' FolderName '\']);
system(['copy VesParam.mat E:\1DWin\' FolderName '\']);
system(['copy NetType.mat E:\1DWin\' FolderName '\']);
system(['copy Structure.mat E:\1DWin\' FolderName '\']);
system(['copy SegName.mat E:\1DWin\' FolderName '\']);
% ����������ӦĿ¼��
system(['copy LinEquDatFile.mat PulseAdapDIR']);
system(['copy ModelParam.mat PulseAdapDIR']);
system(['copy VesParam.mat PulseAdapDIR']);
system(['copy NetType.mat PulseAdapDIR']);
system(['copy Structure.mat PulseAdapDIR']);
system(['copy SegName.mat PulseAdapDIR']);
% ���ڹ�ģ�ϴ�����磬����Ŀ¼���й�����ļ�
if NetTypeID<100
  system('del *.bcs *.in *.his');
end

%% 7. ����run.bat�ű�
runfid=fopen('run.bat','w');
fprintf(runfid,'..\\HemoSim.exe %s.in\r\n@echo off\r\necho \r\n@echo on\r\npause',NetTypeName);
fclose(runfid);
system(['copy run.bat E:\1DWin\' FolderName]);
system(['copy run.bat PulseAdapDIR']);
