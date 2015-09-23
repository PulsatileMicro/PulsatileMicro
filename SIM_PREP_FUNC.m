function [SS_Press,SS_Flow] = SIM_PREP_FUNC(NetTypeID,Boundary,DatMatrix)
global DampInit DampVisc DampE
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP Sparse_SS
global NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE IN_BOUND_TYPE RANDOM_PHASE VEL_PROFILE USE_OPTBOUND BIFURLOSS
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global HUMAN CATMEN HALFSIN IMPULSE INCRAMP DECRAMP SINE MULTI_FREQ_SIN SQUARE HUMAN_ROUTEA_546 HYPER VIT

% 参数意义详见Macro.m
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

% ModelParam=[0 0 1 0 0 0 1 ONED_IMP 0 4 0 0 1e-6 1e-6 1 0 0 1 10 100 1 0]; % 0D血管网络仿真参数
%% Step 1: 准备工作
%%%% 1.1 根据网络类型，调整仿真周期数 %%%%
% TODO: 类型取名
% 仿真nCycle个心动周期
switch NetTypeID
  case {Egg_818_ID,Egg_636_ID,Egg_CAM_ID,Sub_CAM_ID}
    nCycle=12;    % Egg Period=280ms, in total 3360 steps
  otherwise
    nCycle=4;       % Mes Period=800ms, in total 3200 steps
end

%%%% 1.2 Damping因素排除模式: DampFactor %%%%
% DampFactor: DampInit, DampVisc, DampE
DampFactor=DampVisc;
ViscRatio=0.1;
ERatio=0.1;
DampPara=[DampFactor,ViscRatio,ERatio];
DampFactorName=GenDampFactorName(DampFactor,ViscRatio,ERatio);
fprintf('DampType=%s',DampFactorName);

%%%% 1.3 输入信号周期倍数 %%%%
% 周期=周期/nFreqTimes，调整周期
switch NetTypeID
  case {Egg_818_ID,Egg_636_ID,Egg_CAM_ID,Sub_CAM_ID}
    nFreqTimes=1;   % 对于Egg，取其原值
  otherwise
%     nFreqTimes=4;   % 对于Mes，因为所用输入波形采自人体，HR约为75bpm，取其4倍即300bpm作为输入
    nFreqTimes=1;     % 测试
end

% 根据不同血管类型，不同周期倍数，调整信号周期 %%%%
Period=GetPeriod(NetTypeID,ModelParam)/nFreqTimes;
nCycle=nCycle*nFreqTimes;

%%%% 1.4 边界值调节率 %%%%
% 计算边界阻力R时，实际流量大于边界血管处测得的流量，因此需进行调节
switch NetTypeID
  case {Egg_818_ID,Egg_636_ID,Egg_CAM_ID,Sub_CAM_ID}
    OutBoundRatio=1;    % 对于Egg，基本不需要调节，因为其网络出口很接近心脏了
  otherwise
%     OutBoundRatio=100;  % 对于Mes，需要调节，但调节多少无法确定(灵敏度分析TODO)，暂定100
    OutBoundRatio=1;    % 测试
end

%%%% 1.5 仿真步长 %%%%
dt=GetDt(ModelParam(ODESOLVER));

%%%% 1.6 0D方法中调整量纲 %%%%
if ModelParam(ODESOLVER)==RLC_EXP || ModelParam(ODESOLVER)==RLC_IMP ||...
    ModelParam(ODESOLVER)==RC_EXP || ModelParam(ODESOLVER)==RC_IMP
  len_ratio=10;
  mass_ratio=100;
else
  len_ratio=1;
  mass_ratio=1;
end

%%%% 1.7 NetTypeName用于生成.in文件 %%%%
% NetTypeName用于给输出的文件命名
NetTypeName=GetNetTypeName(NetTypeID);

%%%% 1.8 SolverName %%%%
SolverName=GetSolverName(ModelParam(ODESOLVER));

%% Step 2: 稳态血流仿真
% 进行动态仿真前，先仿真稳态血流动力学特性，仿真结果用于为动态模型设置边界条件，主要是输出边界条件
[SS_Press,SS_Flow,SS_Vel,SS_DeltaP,SS_Visc,SS_Hd]=LinEqu(NetTypeID,DatMatrix,Boundary,DampPara);
VesNum=length(SS_Press);

%% Step 3: 分析1D模型边界参数存储矩阵初始化
% 此处边界指的是所有血管段的边界，而非仅仅指血管网络的出入边界
% 所有血管段入口与出口处的边界类型，大小为VesNum*2
BCTypeAll=[];
% 对于分叉、汇聚、接合边界，Bifur为其子血管编号
% 对于输入边界，Bifur为其输入类型编号，如u 0, u 1的0和1
% 对于输出边界，Bifur不起作用
Bifur=zeros(VesNum,4);
% 对于分叉、汇聚、接合边界，BCVal不起作用
% 对于输入边界，BCVal为输入边界的值，如a=PI*1e-4
% 对于输出边界，BCVal为输出边界R或T的值
BCVal=zeros(VesNum,1);
% 边界流速矩阵
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

%%%% 为所有血管设置边界
for i=1:VesNum
  BCType=[];  % 临时向量，存储边界类型
  %%% 分析输入边界
  inInd1=find(From(i)==From); % 查询其他从该段血管起点出发的血管
  inInd2=find(From(i)==To);   % 查询流入该段血管起点的血管
  if length(inInd1)==2
    % 如果有2条血管具有该起点，说明该节点是一个分叉节点
    % 此时必有且仅有一条血管流入该点，即inInd2
    BCType=[BCType 'B'];
    inInd1(inInd1==i)=[]; % 删除本身这条血管的编号
    Bifur(i,1:2)=[inInd2 inInd1];
  elseif length(inInd2)==2
    % 如果有2条血管流入该起点，说明该节点是一个汇聚节点
    BCType=[BCType 'C'];
    Bifur(i,1:2)=inInd2;
  elseif length(inInd2)==0
    % 如果没有血管流入该起点，说明该段血管是一个流入边界
    switch ModelParam(16)
      case 0
        BCType=[BCType 'u'];
        BCVal(i)=SS_Vel(i);  % Velocity in mm/s
        Bifur(i,1:2)=[SS_Hd(i) 2];
      case 1
        BCType=[BCType 'q'];
        BCVal(i)=SS_Vel(i);  % Flow in mm^3/s, 管径在GenBoundInput函数中再乘
        if ModelParam(ODESOLVER)==RLC_EXP || ModelParam(ODESOLVER)==RLC_IMP ||...
            ModelParam(ODESOLVER)==RC_EXP || ModelParam(ODESOLVER)==RC_IMP   % 如果是0D模型
          BCVal(i)=BCVal(i);
        elseif ModelParam(ODESOLVER)==SS || ModelParam(ODESOLVER)==Sparse_SS
          % 如果是稳态模型，边界值设为流量值
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
    % 以上情况均不符，则为接合血管
    BCType=[BCType 'J'];
    Bifur(i,1:2)=inInd2;
  end
  
  %%% 分析输出边界
  outInd1=find(To(i)==From);  % 查询该段血管流向的血管
  outInd2=find(To(i)==To);    % 查询其他到达该段血管终点的血管
  if length(outInd1)==2
    % 如果该段血管流向两根血管，说明该节点是一个分叉节点
    BCType=[BCType 'B'];
    Bifur(i,3:4)=outInd1;
  elseif length(outInd2)==2
    % 如果有2条血管到达该终点，说明该节点是一个汇聚节点
    BCType=[BCType 'C'];
    outInd2(outInd2==i)=[];
    Bifur(i,3:4)=[outInd2 outInd1];
  elseif length(outInd1)==0
    % 如果该终点不是任何血管的起点，说明该节点是一个流出边界
    if ModelParam(ODESOLVER)~=SS && ModelParam(ODESOLVER)~=Sparse_SS % 不是稳态模型，设置R边界
      BCType=[BCType 'R'];
      % 真实的SS_Flow应该更大，所以BCVal应该更小
      % SS_Flow本身不能随OutBoundRatio变化
      BCVal(i)=(SS_Press(i)-SS_DeltaP(i)/2)/SS_Flow(i)/OutBoundRatio*133*60*1e12;
      if ModelParam(ODESOLVER)==RLC_EXP || ModelParam(ODESOLVER)==RLC_IMP ||...
          ModelParam(ODESOLVER)==RC_EXP || ModelParam(ODESOLVER)==RC_IMP
        BCVal(i)=BCVal(i)/1e9*mass_ratio/len_ratio^3;
      end
      % 根据DampFactor调节边界阻力
      % 粘滞度调节时，边界阻力值也必须调节
      % 但是不能在线性方程仿真前调节，否则会导致稳态与动态模型血压水平不一致
      % 是否应该调，还需讨论，可能应该分Visc01和Visc01BD01两类
      if DampFactor==DampVisc
        BCVal(i)=BCVal(i)*ViscRatio;
      end
    else
      % 稳态模型
      ind=find(To(i)==Boundary(:,1));
      if Boundary(ind,2)==0
        % 如果边界类型为压力，则设置血压为边界
        BCType=[BCType 'R'];
        BCVal(i)=SS_Press(i);
      else
        % 如果边界类型为流量，则设置流量为边界，符号为负
        BCType=[BCType 'q'];
        BCVal(i)=-SS_Flow(i)/60/1e12;
      end
    end
  else
    % 以上情况均不符，则为接合血管
    BCType=[BCType 'J'];
    Bifur(i,3:4)=outInd1;
  end
  BCTypeAll=[BCTypeAll;BCType];
end

%% 4. 设置血管参数 %%%%
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
  VesParam(6,:)=4/3;                % alpha=4/3时，等价于0D模型的设置
elseif ModelParam(19)==1
  VesParam(6,:)=1.26;
end
VesParam(7,:)=SS_Visc'*1e-3; % Viscosity
if ModelParam(ODESOLVER)==RC_IMP || ModelParam(ODESOLVER)==RC_EXP ||...
    ModelParam(ODESOLVER)==RLC_IMP || ModelParam(ODESOLVER)==RLC_EXP % 0D模型，调节量纲
  VesParam(1:3,:)=VesParam(1:3,:)*1e3*len_ratio;
  VesParam([4 7],:)=VesParam([4 7],:)*mass_ratio;
end
VesParam(8,:)=SS_Hd;              % Discharge hematocrit
VesParam(9,:)=SegName;            % Segment Name
VesParam(10,:)=From;              % From nodes
VesParam(11,:)=To;                % To nodes
VesParam(12,:)=SS_Vel;
% 设置1D模型的粘弹性血管壁模型参数
if ModelParam(VISCOELAS)==1
  VesParam(13,:)=Eval_Gamma_New(VesParam(2,:),VesParam(5,:));  % Visc part of viscoelasticity
elseif ModelParam(VISCOELAS)==2
  VesParam(14,:)=Eval_Gamma_Old(VesParam(2,:),VesParam(5,:),VesParam(3,:));  % Visc part of viscoelasticity
else
  VesParam(13,:)=0;
  VesParam(14,:)=0;
end
if VesNum>10  % 对于血管网络，1D模型中的q,L参数设为3
  VesParam(15,:)=3;                   % q
  VesParam(16,:)=3;                   % L
else
  VesParam(15,:)=20;                   % q
  VesParam(16,:)=20;                   % L
end
% 去量纲化参数
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
VesParam(22,:)=SegType;             % 血管段类型
VesParam(23,:)=ModelParam(VISCOELAS);

%% 5. 生成边界输入数据文件 %%%%
% 模型从所生成的文件中读取边界波形
% Steady State模型不需要
if ModelParam(ODESOLVER)~=SS && ModelParam(ODESOLVER)~=Sparse_SS
  inFileName=[NetTypeName '_IN'];
  %   outFileName=[NetTypeName '_OUT'];
  outFileName=[NetTypeName '_IN'];
  [i_artinput i_veninput]=...
    GenBoundInput(inFileName,outFileName,BoundData,ModelParam,...
    VesParam,nFreqTimes,0,dt,Period,nCycle,len_ratio,mass_ratio);
end

%     % 确保各种Damping情况下，边界条件一致。主要是出边界R一致
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
%       if BCVal(i)>500 % 大于500的应该是边界阻力值，小于500的应该是边界流量值
%         if ModelParam(ODESOLVER)==RLC_EXP || ModelParam(ODESOLVER)==RLC_IMP ||...
%             ModelParam(ODESOLVER)==RC_EXP || ModelParam(ODESOLVER)==RC_IMP
%           BCVal(i)=BCVal(i)/1e9*mass_ratio/len_ratio^3;
%         end
%       end
%     end

%% 6. 生成1D模型仿真所需读取的.in文件 %%%%
fileName=GenInputFile(NetTypeName,VesParam,BCTypeAll,Bifur,BCVal,ModelParam,DampFactorName,dt,nCycle,Period,len_ratio,mass_ratio);
% 将仿真预处理的信息与仿真结果保存于同一目录，以便Review程序分析仿真结果
save('ModelParam.mat','ModelParam','dt','Period','nCycle','len_ratio','mass_ratio');
save('VesParam.mat','VesParam');
save('NetType.mat','NetTypeName','NetTypeID','Net_546_ID','Net_546_Meas_ID','Egg_818_ID','Net_122_ID','Net_389_ID','Net_913_ID','Egg_CAM_ID','Sub_CAM_ID','Egg_636_ID');
save('Structure.mat','Diam');
save('SegName.mat','SegName');
% 将预处理结果拷贝至指定位置进行仿真
% 根据仿真参数设置目录名，包括1.网络名, 2.求解器类型, 3.Damp因子名, 4.是否Viscoelastic wall 5.心率, 6.边界调节比例
FolderName=[NetTypeName '_' SolverName '_' DampFactorName '_Gamma' int2str(ModelParam(VISCOELAS)) '_bpm' int2str(60/Period)...
  '_BDFlow' int2str(OutBoundRatio) '_BifurLoss' int2str(ModelParam(BIFURLOSS))];
system(['mkdir E:\1DWin\' FolderName]);
system(['mkdir PulseAdapDIR']);
system(['copy ' NetTypeName '*.* E:\Projects\HemoSim\HemoSim\']);
system(['copy ' NetTypeName '*.* E:\Projects\HemoSim_Sundials2.6\HemoSim\']);
system(['copy ' NetTypeName '*.* E:\1DWin\' FolderName '\']);
system(['copy ' NetTypeName '*.* PulseAdapDIR\']);
% 拷贝到运行目录下
system(['copy LinEquDatFile.mat E:\1DWin\' FolderName '\']);
system(['copy ModelParam.mat E:\1DWin\' FolderName '\']);
system(['copy VesParam.mat E:\1DWin\' FolderName '\']);
system(['copy NetType.mat E:\1DWin\' FolderName '\']);
system(['copy Structure.mat E:\1DWin\' FolderName '\']);
system(['copy SegName.mat E:\1DWin\' FolderName '\']);
% 拷贝到自适应目录下
system(['copy LinEquDatFile.mat PulseAdapDIR']);
system(['copy ModelParam.mat PulseAdapDIR']);
system(['copy VesParam.mat PulseAdapDIR']);
system(['copy NetType.mat PulseAdapDIR']);
system(['copy Structure.mat PulseAdapDIR']);
system(['copy SegName.mat PulseAdapDIR']);
% 对于规模较大的网络，避免目录中有过多的文件
if NetTypeID<100
  system('del *.bcs *.in *.his');
end

%% 7. 生成run.bat脚本
runfid=fopen('run.bat','w');
fprintf(runfid,'..\\HemoSim.exe %s.in\r\n@echo off\r\necho \r\n@echo on\r\npause',NetTypeName);
fclose(runfid);
system(['copy run.bat E:\1DWin\' FolderName]);
system(['copy run.bat PulseAdapDIR']);
