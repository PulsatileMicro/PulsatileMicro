function [InvSeg,AllBoundFlow] = BoundOptSS(DampFactor,NetTypeID,VesNum,SegType,From,To,Diam,Len,Visc,Boundary)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN

% 边界条件
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% PSO算法参数设置
nItNum=20;                    % 迭代次数
nPSO=25;                      % 粒子群中粒子个数
nParDim=length(BoundFlow);    % 每个粒子的维度. 此处设置为length(BoundFlow),但实际优化时会去掉主输入和主输出
AllBoundFlow=zeros(nParDim,nPSO,nItNum);  % 粒子位置xi
v=zeros(nParDim,nPSO,nItNum);             % 粒子速度vi
AllBoundFlowLimit=zeros(nParDim,2);       % 粒子中，每个元素的最大值最小值
InvSeg=VesNum*ones(nPSO,nItNum);  % 适应度结果——逆流血管段数，初值为VesNum段
InvSegVel=1e100*ones(nPSO,nItNum);% 备用适应度结果——逆流血管的流速和
Pg=zeros(nParDim,1);              % 某次搜索全局最优
Pi=zeros(nParDim,nPSO);           % 各粒子历史最优
AllMeanP=zeros(VesNum,nItNum);
AllMeanFlow=zeros(VesNum,nItNum);

% 边界血管方向修正 %%%%
% 依据边界条件的符号修正
% 血流为正数->入边界->边界节点在From类中
% 血流为负数->出边界->边界节点在To类中
[From,To]=OriNodeModify(Boundary,From,To);
Qref=0.001; % TODO: Qref没用，这个程序需要重构
% Boundary(:,4)=0.45;
BHd=BoundaryInput(Boundary,From,To,Qref);

% 单位调整
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s

% 设置主输入下标(在边界向量中的下标)，主输入边界不参与优化
switch NetTypeID
  case {Net_546_Meas_ID,Net_546_ID}
    MainBoundInd=[23 28];
  case Egg_818_ID
    MainBoundInd=[11 20 27];
  case Net_389_ID
    MainBoundInd=[6 7];
  case Net_913_ID
    MainBoundInd=[1 34];
  case Egg_636_ID
    MainBoundInd=[1 13];
end
% 初始化粒子群
% 每个粒子中的每个元素初始化为随机数
for i=1:nPSO
  for j=1:length(BoundFlow)
    ind=find(j==MainBoundInd);
    if ~isempty(ind)
      AllBoundFlow(MainBoundInd(ind),i,1)=BoundFlow(MainBoundInd(ind));
    else
      % 对于非主输入、主输出边界，设置随机初始值
      % TODO: 如何设？
      AllBoundFlow(j,i,1)=BoundFlow(j)*rand*1e2;
    end
  end
end
% 设置粒子中元素的上下限，边界的血流方向保持不变
AllBoundFlowLimit(:,1)=mean(AllBoundFlow(:,:,1),2)*1e-3;
AllBoundFlowLimit(:,2)=mean(AllBoundFlow(:,:,1),2)*1e3;
AllBoundFlowLimit(MainBoundInd,1)=AllBoundFlow(MainBoundInd,1,1);
AllBoundFlowLimit(MainBoundInd,2)=AllBoundFlow(MainBoundInd,1,1);
% 如果边界为负，则元素上下限对调，确保比较大小时逻辑正确
for i=1:length(BoundFlow)
  if BoundFlow(i)<0
    temp=AllBoundFlowLimit(i,1);
    AllBoundFlowLimit(i,1)=AllBoundFlowLimit(i,2);
    AllBoundFlowLimit(i,2)=temp;
  end
end
% 速度、位置更新公式参数
c1=4; % 认知
c2=2; % 社会
r=1;  % 约束因子
w=0.729;  % 惯性权重
CTime=5;
for m=1:nItNum  % 迭代nItNum次
  m
  c1=c1*0.8;
  c2=c2*1.2;
  for k=1:nPSO  % 对于每一个粒子
    T=0;
    Visc=2*ones(VesNum,1)*1e-3;       % 初始粘滞度
    Porder=1:VesNum;Norder=VesNum:1;  % 初始Hd计算顺序
    errorFlag=0;
    while T<CTime   % 循环1，Visc反馈
      T=T+1
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
      
      if NetTypeID~=Net_546_ID && NetTypeID~=Net_122_ID
        %%粘滞度反馈
        umDiam=Diam.*1e6;   %um
        for i=1:VesNum
          Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
          if DampFactor==Visc01
            Visc(i)=Visc(i)/10;
          end
        end
        if T>1
          Visc=0.5*(Visc+DebugVisc(:,T-1));
        end
      end
      
      % 记录每次迭代的仿真结果
      DebugVisc(:,T)=Visc;
      DebugHd(:,T)=Hd;
      DebugPressure(:,T)=MeanP;
      DebugFlow(:,T)=MeanFlowNew;
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
end