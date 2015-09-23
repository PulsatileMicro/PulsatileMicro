%% 优化血管网络边界条件
function [AllBoundFlow,InvSeg] = BOUND_OPT_FUNC(NetTypeID,Boundary,DatMatrix)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global NOT_OPT_PARA OPT_PARA OPT_BOUND_FLOW OPT_BOUND_META
OptType=OPT_BOUND_FLOW;
%% 1. 根据不同类型的血管网络确定主输入、主输出边界
% 所确定的编号为其序号，而非血管段名称
switch NetTypeID
  case {Net_546_ID,Net_546_Meas_ID}
    MainIn=1;
    MainOut=330;
  case Net_389_ID
    MainIn=100;
    MainOut=87;
  case Net_913_ID
    MainIn=376;
    MainOut=395;
  case Net_122_ID
    % TODO
  case Egg_818_ID
    MainIn=[761,765,216,710];
    MainOut=357;
  case Egg_636_ID
    MainIn=166;
    MainOut=2;
  case Egg_CAM_ID
    % TODO
  case Sub_CAM_ID
    % TODO
end

%% 2. 数据解包
% 血管段数据
SegName=DatMatrix(:,1);
From=DatMatrix(:,2);
To=DatMatrix(:,3);
Len=DatMatrix(:,4);
Diam=DatMatrix(:,5);
WallTh=DatMatrix(:,6);
SegType=DatMatrix(:,7);
Visc=DatMatrix(:,8);
E=DatMatrix(:,9);
VesNum=length(Len);

% 边界数据
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
% BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);
BoundMainInOut=zeros(length(BoundNode),1);

for i=1:length(BoundNode)
  % 查找边界节点为入口节点的血管的编号，为了找到主输入
  SegInd=find(From==BoundNode(i));
  if BoundType(i)==0
    % 压力边界
    BoundFlow(i)=Boundary(i,3);
    BoundMainInOut(i)=1;
  elseif ~isempty(SegInd)
    % 表明该节点为输入节点
    if ~isempty(find(MainIn==SegInd,1))
      % 查到主输入
      BoundFlow(i)=Boundary(i,3);
      BoundMainInOut(i)=1;
    else
      % 一般输入
%       BoundFlow(i)=(0.4*Diam(SegInd)-1.9)*0.25*pi.*Diam(SegInd).^2/1e3*60;
      BoundFlow(i)=30;
    end
  else
    SegInd=find(To==BoundNode(i));
    %  输出节点
%     BoundFlow(i)=-(0.4*Diam(SegInd)-1.9)*0.25*pi.*Diam(SegInd).^2/1e3*60;
    BoundFlow(i)=-30;  % nL/min
  end
end

% 单位调整
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s

% 边界血管方向修正 %%%%
% 依据边界条件的符号修正
% 血流为正数->入边界->边界节点在From类中
% 血流为负数->出边界->边界节点在To类中
[From,To]=OriNodeModify(Boundary,From,To);

Qref=0.001; % TODO: Qref没用，这个程序需要重构
% Boundary(:,4)=0.45;
BHd=BoundaryInput(Boundary,From,To,Qref);

%% 粘滞度循环
% 循环计数
Loop1_Cnt=0;
Loop1_Visc_Limit=20;

DebugVisc=zeros(VesNum,Loop1_Visc_Limit);
DebugHd=zeros(VesNum,Loop1_Visc_Limit);
DebugPressure=zeros(VesNum,Loop1_Visc_Limit);
DebugFlow=zeros(VesNum,Loop1_Visc_Limit);
Porder=1:VesNum;Norder=VesNum:1;  % 初始Hd计算顺序

%% 粒子群算法参数
nItNum=20;                    % 迭代次数
nPSO=25;                      % 粒子群中粒子个数
nParDim=length(BoundFlow);    % 每个粒子的维度
AllBoundFlow=zeros(nParDim,nPSO,nItNum);  % 所有粒子xi
AllBoundFlowLimit=zeros(nParDim,2);       % 粒子的最大值最小值
v=zeros(nParDim,nPSO,nItNum);             % 粒子速度vi

for i=1:nPSO
  % 初始化粒子群
  for j=1:length(BoundFlow)
    if BoundMainInOut(j)
      AllBoundFlow(j,i,1)=BoundFlow(j);
    else
      AllBoundFlow(j,i,1)=BoundFlow(j)*rand;
    end
  end
end
% 设置粒子中元素的上下限
% 边界的血流方向保持不变
AllBoundFlowLimit(:,1)=AllBoundFlow(:,1,1)*1e-3;
AllBoundFlowLimit(:,2)=AllBoundFlow(:,1,1)*1e3;
MainInOutIndex=find(BoundMainInOut);
for i=1:length(MainInOutIndex)
  AllBoundFlowLimit(MainInOutIndex(i),:)=AllBoundFlow(MainInOutIndex(i),1,1);
end
for i=1:length(BoundFlow)
  if BoundFlow(i)<0
    temp=AllBoundFlowLimit(i,1);
    AllBoundFlowLimit(i,1)=AllBoundFlowLimit(i,2);
    AllBoundFlowLimit(i,2)=temp;
  end
end
c1=2; % 认知
c2=2; % 社会
r=1;  % 约束因子
w=0.729;  % 惯性权重
InvSeg=VesNum*ones(nPSO,nItNum);  % 适应度结果——逆流血管段数，初值为VesNum段
InvSegVel=1e100*ones(nPSO,nItNum);% 备用适应度结果——逆流血管的流速和
Pg=zeros(nParDim,1);              % 某次搜索全局最优
Pi=zeros(nParDim,nPSO);           % 各粒子历史最优
AllMeanP=zeros(VesNum,nItNum);
AllMeanFlow=zeros(VesNum,nItNum);

for m=1:nItNum  % 迭代nItNum次
%   c1=c1*0.8;
%   c2=c2*1.2;
  for k=1:nPSO  % 对于每一个粒子
    fprintf('Iter=%d, PSOid=%d ViscLoopCnt=%d\n',m,k,Loop1_Cnt);
    Loop1_Cnt=0;
    Loop1_Visc_Num=Loop1_Visc_Limit;
    Visc=2*ones(VesNum,1)*1e-3;
    Porder=1:VesNum;Norder=VesNum:1;  % 初始Hd计算顺序
    errorFlag=0;
    while Loop1_Cnt<Loop1_Visc_Num   % 循环1，Visc反馈
      Loop1_Cnt=Loop1_Cnt+1;
      % 线性方程求解模块
      [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,AllBoundFlow(:,k,m),VesNum,Diam,Visc,Len,From,To);
      % 逆流处理模块
      [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju);
      %     FromNew=From;ToNew=To;MeanFlowNew=MeanFlow;DeltaPNew=DeltaP;
      % 调整Hd计算顺序
      [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,2,Eju);
      if Eju==1
        errorFlag=1;
        break;
      end
      % 计算Hd
      [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam,MeanFlowNew,Eju);
      % 粘滞度反馈
      umDiam=Diam.*1e6;   % um
      for i=1:VesNum
        Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
      end
      DebugVisc(:,Loop1_Cnt)=Visc;
      % 循环进入第二次后(Loop1_Cnt>1)，开始校正Visc
      % 采用SOR方法校正，避免出现粘滞度振荡的情况
      ViscModType=1;  % 0：不修正；1：对Bifurcation修正；2：对All segments修正
      Alpha=0.5;
      DebugVisc=modifyViscosity(Porder,FromNew,ToNew,DebugVisc,Loop1_Cnt,ViscModType,Alpha,Eju);
      Visc=DebugVisc(:,Loop1_Cnt);
      
      %%%% 判断循环结束模块 %%%%
      % LoopOutType:0-正常；1-收敛；2-达到步长未收敛；3-visc错误；4-线性方程错误
      [ViscMAE,Loop1OutType,Loop1_Visc_Num]=ViscLoopTerminator(DebugVisc*1e3,Loop1_Visc_Num,Loop1_Cnt,Eju);
    end
    if errorFlag~=1
      InvSeg(k,m)=length(find(MeanFlow<0));   % 每一个粒子对应的适应度函数的结果
      MeanVel=MeanFlow./(0.25*pi*Diam.^2);
      InvSegVel(k,m)=-sum(MeanVel(MeanVel<0));
    end
  end
  % 优化目标为逆流血管段数
  [Value1,Ind1]=min(InvSeg);          % Ind为全部粒子搜索到的历史最优值Pg的序号
  [Value2,Ind2]=min(Value1);
  Pg=AllBoundFlow(:,Ind1(Ind2),Ind2);         % Pg为本次(m)搜索中的全局最优值（是否应该取历史搜索中的全局最优值？）
  % 寻找历史中每个粒子的最优值
  for k=1:nPSO
    [Value,Ind]=min(InvSeg(k,:));   % Ind为第k个粒子的历史最优值的序号
    Pi(:,k)=AllBoundFlow(:,k,Ind);  % Pi(:,k)为第k个粒子的历史最优值
    v(:,k,m+1)=w.*v(:,k,m)+c1*rand*(Pi(:,k)-AllBoundFlow(:,k,m))+c2*rand*(Pg-AllBoundFlow(:,k,m));
    AllBoundFlow(:,k,m+1)=AllBoundFlow(:,k,m)+r*v(:,k,m+1);
    for i=1:nParDim
      if AllBoundFlow(i,k,m+1)<AllBoundFlowLimit(i,1)
        AllBoundFlow(i,k,m+1)=AllBoundFlowLimit(i,1);
      elseif AllBoundFlow(i,k,m+1)>AllBoundFlowLimit(i,2)
        AllBoundFlow(i,k,m+1)=AllBoundFlowLimit(i,2);
      end
    end
  end
  % 记录本次迭代中的最小MeanP和MeanFlow
  AllMeanP(:,m)=MeanP;
  AllMeanFlow(:,m)=MeanFlow;
end

Vel=MeanFlowNew/1e12/60./(pi*0.25*Diam.^2)*1e3;