function [AllFlowX,AllPressX,NormAllFlowX,NormAllPressX,CutFreqFlow,CutFreqPress,AllDeltaPX,AllRatioPX]=...
  LinEquFreq(VesNum,From,To,Diam,Len,Visc,WallTh,E,Boundary,SS_Press,SS_Flow,SS_DeltaP,DampFactor,ViscRatio,InputType,f,df)
%% 血管网络的频域分析
%%%% 边界数据
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% 单位调整
Len=Len*1e-6;       %m
Diam=Diam*1e-6;     %m
Visc=Visc*1e-3;     %Pa.s
WallTh=WallTh*1e-6; %m

% 去量纲，保证线性方程求解效率(TODO: 未实现)
len_ratio=1;
mass_ratio=1;
Len=Len*len_ratio;
Diam=Diam*len_ratio;
Visc=Visc*mass_ratio;
WallTh=WallTh*len_ratio;
E=E*mass_ratio;

% 计算阻力和顺应性
R=(8.*Visc.*Len)./(pi*(Diam./2).^4);
C=3*pi*(Diam./2).^3.*Len./(2.*E.*WallTh);

DeltaPX=zeros(VesNum,1);PressXX=DeltaPX;RatioPX=DeltaPX;RatioP2=DeltaPX;
% 准备求解线性方程组
% AllNodes=union(From,To);
% JNum=length(AllNodes);
% JMat=zeros(JNum+VesNum,JNum+VesNum);  % 左边矩阵
% RHS=zeros(JNum+VesNum,1); % 右边矩阵，顺序为Q1->QVesNum,P1->PJNum
AllNodes=union(From,To);
JNum=length(AllNodes);
% f(1)=0.001;
% f=[0 5 10];
AllPressX=zeros(length(f),VesNum);AllFlowX=AllPressX;AllDeltaPX=AllPressX;AllRatioPX=AllPressX;
for k=1:length(f)
  fprintf(1,'Freq=%.3fHz\n',f(k));
  s=1i*2*pi*f(k);
  X=LinEquFreqSolver(s,R,C,VesNum,JNum,AllNodes,From,To,SS_Press,SS_Flow,SS_DeltaP,len_ratio,mass_ratio,DampFactor,ViscRatio,InputType,Boundary);
  absX=abs(X);
  FlowX=absX(1:VesNum);
  PressX=absX(VesNum+1:end);
  
  % 计算每段血管的平均血压和首尾压差
  for j=1:VesNum
    InIndex=find(From(j)==AllNodes);
    OutIndex=find(To(j)==AllNodes);
    %     if PressX(InIndex)>=PressX(OutIndex)
    DeltaPX(j)=PressX(InIndex)-PressX(OutIndex);
    RatioPX(j)=PressX(OutIndex)./PressX(InIndex);
    PressXX(j)=PressX(InIndex)-DeltaPX(j)/2;
    %     else
    %       DeltaPX(j)=PressX(OutIndex)-PressX(InIndex);
    %       RatioPX(j)=PressX(InIndex)./PressX(OutIndex);
    %       PressXX(j)=PressX(OutIndex)-DeltaPX(j)/2;
    %     end
  end
  AllFlowX(k,:)=FlowX;
  AllPressX(k,:)=PressXX;
  AllDeltaPX(k,:)=DeltaPX;
  AllRatioPX(k,:)=RatioPX;
end

% 计算每段血管的截止频率
CutFreqFlow=zeros(VesNum,1);CutFreqPress=zeros(VesNum,1);
% 每个频率分量上的值都相对入口血管做归一化
% TODO: 入口血管不一定是编号1
for i=1:length(f)
  AllPressX(i,:)=AllPressX(i,:)./AllPressX(i,1);
  AllFlowX(i,:)=AllFlowX(i,:)./AllFlowX(i,1);
end
for i=1:VesNum
  % Norm的意义：以0Hz为参照，随着频率的增高，衰减如何变化
  NormAllFlowX(:,i)=AllFlowX(:,i)./(AllFlowX(1,i));
  NormAllPressX(:,i)=AllPressX(:,i)./(AllPressX(1,i));
  NormAllDeltaPX(:,i)=abs(AllDeltaPX(:,i))./(abs(AllDeltaPX(1,i)));
  NormAllRatioPX(:,i)=abs(AllRatioPX(:,i))./(abs(AllRatioPX(1,i)));
  
  [Value Ind]=min(abs(NormAllFlowX(:,i)-0.707));
  if Ind~=1
    CutFreqFlow(i)=(Ind-1)*df;
  else
    CutFreqFlow(i)=Ind*df;
  end
  [Value Ind]=min(abs(NormAllPressX(:,i)-0.707));
  if Ind~=1
    CutFreqPress(i)=(Ind-1)*df;
  else
    CutFreqPress(i)=Ind*df;
  end
end

