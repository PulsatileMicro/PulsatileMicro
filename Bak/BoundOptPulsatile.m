function [InvSeg,DiffSeg,AllBoundFlow] = BoundOptPulsatile(NetTypeID,SegType,VesNum,From,To,Diam,Len,Visc,WallTh,E,Boundary,SS_Press,SS_Flow,SS_DeltaP)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit
%%%% 边界数据
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% 设置主输入下标，主输入边界不参与优化
switch NetTypeID
  case {Net_546_Meas_ID,Net_546_ID}
    MainBoundInd=[23 28];
  case Egg_818_ID
    %         MainInputInd
end

% 频率分析参数
df=2;    % 步长
f=0:df:10;

% PSO算法参数设置
nItNum=7;                                   % 迭代次数
nPSO=10;                                    % 粒子群中粒子个数
nParDim=length(BoundFlow);                  % 每个粒子的维度. 此处设置为length(BoundFlow),但实际优化时会去掉主输入和主输出
AllBoundFlow=zeros(nParDim,nPSO,nItNum);    % 粒子位置xi(即粒子的值)
v=zeros(nParDim,nPSO,nItNum);               % 粒子速度vi
AllBoundFlowLimit=zeros(nParDim,2);         % 粒子中，每个元素的最大值最小值，防止粒子随机变化过大，或改变方向
InvSeg=zeros(VesNum,length(f),nPSO,nItNum); % 每次迭代，每个频率分量，粒子群中每个粒子作用下，逆流血管编号
DiffSeg=1e100*ones(nPSO,nItNum);
Pg=zeros(nParDim,1);                        % 某次搜索全局最优
Pi=zeros(nParDim,nPSO);                     % 各粒子历次、历史最优

% 初始化粒子群
% 第一次迭代所需的粒子群中，每个粒子的每个元素初始化为随机数
for i=1:nPSO
  for j=1:length(BoundFlow)
    ind=find(j==MainBoundInd);
    if ~isempty(ind)
      AllBoundFlow(MainBoundInd(ind),i,1)=BoundFlow(MainBoundInd(ind));
    else
      % 对于非主输入、主输出边界，设置随机初始值
      % TODO: 如何设？目前设为原边界值*随机数*10
      AllBoundFlow(j,i,1)=BoundFlow(j)*rand;
    end
  end
end
% 设置粒子中元素的上下限，边界的血流方向保持不变
AllBoundFlowLimit(:,1)=mean(AllBoundFlow(:,:,1),2)*1e-3;
AllBoundFlowLimit(:,2)=mean(AllBoundFlow(:,:,1),2)*5;
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

% 单位调整
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s
WallTh=WallTh*1e-6; %m

% 去量纲，保证线性方程求解效率
% m_ratio=1e6;
% kg_ratio=1e2;
m_ratio=1;
kg_ratio=1;
Len=Len*m_ratio;
Diam=Diam*m_ratio;
Visc=Visc*kg_ratio;
WallTh=WallTh*m_ratio;
E=E*kg_ratio;

% 准备求解线性方程组
DeltaPX=zeros(VesNum,1);PressXX=DeltaPX;RatioPX=DeltaPX;RatioP2=DeltaPX;
% JMat=zeros(JNum+VesNum,JNum+VesNum);  % 左边矩阵
% RHS=zeros(JNum+VesNum,1); % 右边矩阵，顺序为Q1->QVesNum,P1->PJNum
AllNodes=union(From,To);
JNum=length(AllNodes);

% 速度、位置更新公式参数
c1=2; % 认知
c2=2; % 社会
r=1;  % 约束因子
w=0.729;  % 惯性权重

for m=1:nItNum
  fprintf(1,'ItNum=%d\n',m);
%   c1=c1*0.8;
%   c2=c2*1.2;
  for i=1:nPSO
    fprintf(1,'  nPSO=%d\n',i);
    % 以AllBoundFlow为边界，仿真稳态血流动力学特性，从而得到基于随机边界的
    % SS_Press,SS_Flow和SS_DeltaP
    Boundary(:,3)=AllBoundFlow(:,i,m);
    [SS_Press,SS_Flow,SS_Vel,SS_DeltaP,SS_Visc,SS_Hd]=...
      LinEqu(NetTypeID,DampInit,VesNum,SegType,From,To,Diam,Len,Visc,Boundary);
    % 更新阻力和顺应性参数
    R=(8.*Visc.*Len)./(pi*(Diam./2).^4);
    C=3*pi*(Diam./2).^3.*Len./(2.*E.*WallTh);
    AllPressX=zeros(length(f),VesNum);AllFlowX=AllPressX;AllDeltaPX=AllPressX;AllRatioPX=AllPressX;
    for k=1:length(f)
      fprintf(1,'    Freq=%.3fHz\n',f(k));
      s=1i*2*pi*f(k);
      % 注：LinEquFreqSolver函数虽然传入了SS_Press,SS_Flow和SS_DeltaP，但实际上只用到了这几个向量中的边界血管值
      X=LinEquFreqSolver(s,R,C,VesNum,JNum,AllNodes,From,To,SS_Press,SS_Flow,SS_DeltaP,m_ratio,kg_ratio);
      absX=abs(X);  % 取模
      FlowX=absX(1:VesNum);
      PressX=absX(VesNum+1:end);
      
      % 计算每段血管的平均血压和首尾压差
      for j=1:VesNum
        InIndex=find(From(j)==AllNodes);
        OutIndex=find(To(j)==AllNodes);
%         if PressX(InIndex)>=PressX(OutIndex)
          DeltaPX(j)=PressX(InIndex)-PressX(OutIndex);
          RatioPX(j)=PressX(OutIndex)./PressX(InIndex);
          PressXX(j)=PressX(InIndex)-DeltaPX(j)/2;
%         else
%           DeltaPX(j)=PressX(OutIndex)-PressX(InIndex);
%           RatioPX(j)=PressX(InIndex)./PressX(OutIndex);
%           PressXX(j)=PressX(OutIndex)-DeltaPX(j)/2;
%         end
      end
      AllFlowX(k,:)=FlowX;
      AllPressX(k,:)=PressXX;
      AllDeltaPX(k,:)=DeltaPX;
      AllRatioPX(k,:)=RatioPX;
      % 规定：逆流为1，正流为0
      InvInd=find(AllDeltaPX(k,:)<0); % 找到第k个频率的仿真结果中，血流方向为负的血管段编号
      InvSeg(InvInd,k,i,m)=1;         % 记录第k个频率的仿真结果中逆流的血管
    end
    % 所有频率跑完后，要把结果保存成一个值，DiffSeg(nPSO,nItNum)
%     DiffSeg(i,m)=mean(std(InvSeg(:,:,i,m),0,2));  % 考虑每个频率下，血流方向的变异性
    DiffSeg(i,m)=sum(sum(InvSeg(:,:,i,m),2));
  end
  % 优化目标为逆流血管段数
  [Value1,Ind1]=min(DiffSeg);          % Ind为全部粒子搜索到的历史最优值Pg的序号
  [Value2,Ind2]=min(Value1);
  Pg=AllBoundFlow(:,Ind1(Ind2),Ind2);    % Pg为本次(m)搜索中的全局最优值（是否应该取历史搜索中的全局最优值？）
  % 寻找历史中每个粒子的最优值
  for k=1:nPSO
    [Value,Ind]=min(DiffSeg(k,:));   % Ind为第k个粒子的历史最优值的序号
    Pi(:,k)=AllBoundFlow(:,k,Ind);  % Pi(:,k)为第k个粒子的历史最优值
    v(:,k,m+1)=w.*v(:,k,m)+c1*rand*(Pi(:,k)-AllBoundFlow(:,k,m))+c2*rand*(Pg-AllBoundFlow(:,k,m));
    AllBoundFlow(:,k,m+1)=AllBoundFlow(:,k,m)+r*v(:,k,m+1);
    % 限制边界的随机变动范围
    for i=1:nParDim
      if AllBoundFlow(i,k,m+1)<AllBoundFlowLimit(i,1)
        AllBoundFlow(i,k,m+1)=AllBoundFlowLimit(i,1);
      elseif AllBoundFlow(i,k,m+1)>AllBoundFlowLimit(i,2)
        AllBoundFlow(i,k,m+1)=AllBoundFlowLimit(i,2);
      end
    end
  end
end
end
