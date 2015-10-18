%% 仿真55段人体动脉血管的程序
clear;clc;close all;
Macro;
% 运行频域分析
RunFreqAnal=1;
% 1 - dt
% 2 - Period
% 3 - PeriodNum
% 4 - ODESolver 1:ONED_IMP 8:RC_IMP
ModelParam=[1e-3 0.8 8 1]; % 血管网络仿真参数
NetTypeName='Systemic';
% 读取人体血管结构与拓扑数据
[num,txt,raw]=xlsread('SysData.xlsx','A2:S56');
SegName=num(:,1);
Len=num(:,3);
% 第4列为well-matched diam. 第5列为原始diam.
% Diam=2*num(:,4);
Diam=2*num(:,5);
WallTh=num(:,6);
E=num(:,7);
OutR=num(:,8);
OutC=num(:,10);
ConMat=num(:,11:14);
From=num(:,17);
To=num(:,18);
Order=num(:,19);
% 第14列为well-matched diam. 第15列为原始diam.
% PWV=num(:,15);
PWV=num(:,16);
Eh=PWV.^2.*Diam*1e-2*1050;
PWV_E=1e5;
PWV_h=Eh/PWV_E;
VesNum=length(SegName);

%% 血管参数改变
Diam=Diam;

%% 生成仿真文件
% 1D模型边界参数存储矩阵初始化
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

%% 为所有血管设置边界
for i=1:VesNum
  BCType=[];  % 临时向量，存储边界类型
  if ConMat(i,1)==-1  % 输入边界
    BCType=[BCType 'q'];
    Bifur(i,1:2)=2;
  else
    BCType=[BCType 'B'];  % 对于整体动脉循环而言，只有B，没有C
    Bifur(i,1:2)=ConMat(i,1:2);
  end
  
  if ConMat(i,3)==-2  % 输出边界
%     BCType=[BCType 'R'];
%     Bifur(i,3:4)=OutR(i)*1e10;
    BCType=[BCType 'W'];
    Bifur(i,3:4)=[OutC(i)*1e-10 OutR(i)*1e10];
  else
    BCType=[BCType 'B'];
    Bifur(i,3:4)=ConMat(i,3:4);
  end
  BCTypeAll=[BCTypeAll;BCType];
end

%% %% 设置血管参数 %%%%
VesParam=zeros(22,VesNum);
VesParam(1,:)=Len'*1e-2;            % Length
VesParam(2,:)=Diam'*1e-2;           % Diameter
% VesParam(3,:)=WallTh'*1e-2;       % Wall thickness
% VesParam(4,:)=E*1e6;
VesParam(3,:)=PWV_h;                % Wall thickness
VesParam(4,:)=PWV_E;
VesParam(5,:)=1e4;                  % Viscous wall modulus
VesParam(6,:)=1.1;                  % Alpha
VesParam(7,:)=4e-3;                 % Viscosity
VesParam(8,:)=0.45;                 % Discharge hematocrit
VesParam(9,:)=SegName;              % Segment Name
VesParam(10,:)=[];                  % From nodes
VesParam(11,:)=[];                  % To nodes
VesParam(12,:)=[];
% 设置1D模型的粘弹性血管壁模型参数
VesParam(13,:)=0;
VesParam(14,:)=0;
VesParam(15,:)=8;                   % q
VesParam(16,:)=8;                   % L
% Scale Parameters
VesParam(17,:)=1;
VesParam(18,:)=1;
VesParam(19,:)=1;
VesParam(20,:)=3;                   % Number of history points
VesParam(21,:)=0;                   % Taper Rate

%% 统计R,C,L的值
R=8*VesParam(7,:).*VesParam(1,:)./(pi.*(VesParam(2,:)/2).^4);
C=3*pi*(VesParam(2,:)/2).^3.*VesParam(1,:)./(2.*Eh');

%% 生成边界输入数据文件 %%%%
% 模型从所生成的文件中读取边界波形
inFileName=[NetTypeName '_IN'];
i_artinput=GenBoundInput_Systemic(inFileName,ModelParam);

fileName=GenInputFile_Systemic(NetTypeName,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);
system(['mkdir C:\1DWin\' NetTypeName]);
system(['copy ' NetTypeName '*.* C:\1DWin\' NetTypeName '\']);

if RunFreqAnal==1
  AllNodes=union(From,To);
  JNum=length(AllNodes);
  df=0.5;    % 频率步长
  f=0:df:10;
  f(1)=0.01;
  AllPressX=zeros(length(f),VesNum);AllFlowX=AllPressX;AllDeltaPX=AllPressX;AllRatioPX=AllPressX;
  for k=1:length(f)
    fprintf(1,'Freq=%.3fHz\n',f(k));
    s=1i*2*pi*f(k);
    X=LinEquFreqSolver_systemic(s,R,C,VesNum,JNum,AllNodes,From,To,OutR,OutC);
    absX=abs(X);
    FlowX=absX(1:VesNum);
    PressX=absX(VesNum+1:end);
    
    % 计算每段血管的平均血压和首尾压差
    for j=1:VesNum
      InIndex=find(From(j)==AllNodes);
      OutIndex=find(To(j)==AllNodes);
      DeltaPX(j)=PressX(InIndex)-PressX(OutIndex);
      RatioPX(j)=PressX(OutIndex)./PressX(InIndex);
      PressXX(j)=PressX(InIndex)-DeltaPX(j)/2;
    end
    AllFlowX(k,:)=FlowX;
    AllPressX(k,:)=PressXX;
    AllDeltaPX(k,:)=DeltaPX;
    AllRatioPX(k,:)=RatioPX;
  end
  
  RawPressX=AllPressX;
  RawFlowX=AllFlowX;
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
  end
end

MeanNormPressX=zeros(length(f),max(Order));
MeanNormFlowX=MeanNormPressX;
StdNormPressX=MeanNormFlowX;
StdNormFlowX=MeanNormFlowX;
AllOrder=1:max(Order);
for j=1:length(f)
  for i=1:length(AllOrder)
    ind=find((Order==AllOrder(i)));
    MeanNormPressX(j,i)=mean(NormAllPressX(j,ind));
    StdNormPressX(j,i)=std(NormAllPressX(j,ind));
%     DelInd=find(NormAllFlowX(j,ind)>1 | NormAllFlowX(j,ind)<0);
%     if ~isempty(DelInd)
%       ind(DelInd)=[];
%     end
    MeanNormFlowX(j,i)=mean(NormAllFlowX(j,ind));
    StdNormFlowX(j,i)=std(NormAllFlowX(j,ind));
  end
end
MeanNormPressX=MeanNormPressX';
StdNormPressX=StdNormPressX';
MeanNormFlowX=MeanNormFlowX';
StdNormFlowX=StdNormFlowX';
