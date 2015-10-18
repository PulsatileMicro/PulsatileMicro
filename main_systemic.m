%% ����55�����嶯��Ѫ�ܵĳ���
clear;clc;close all;
Macro;
% ����Ƶ�����
RunFreqAnal=1;
% 1 - dt
% 2 - Period
% 3 - PeriodNum
% 4 - ODESolver 1:ONED_IMP 8:RC_IMP
ModelParam=[1e-3 0.8 8 1]; % Ѫ������������
NetTypeName='Systemic';
% ��ȡ����Ѫ�ܽṹ����������
[num,txt,raw]=xlsread('SysData.xlsx','A2:S56');
SegName=num(:,1);
Len=num(:,3);
% ��4��Ϊwell-matched diam. ��5��Ϊԭʼdiam.
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
% ��14��Ϊwell-matched diam. ��15��Ϊԭʼdiam.
% PWV=num(:,15);
PWV=num(:,16);
Eh=PWV.^2.*Diam*1e-2*1050;
PWV_E=1e5;
PWV_h=Eh/PWV_E;
VesNum=length(SegName);

%% Ѫ�ܲ����ı�
Diam=Diam;

%% ���ɷ����ļ�
% 1Dģ�ͱ߽�����洢�����ʼ��
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

%% Ϊ����Ѫ�����ñ߽�
for i=1:VesNum
  BCType=[];  % ��ʱ�������洢�߽�����
  if ConMat(i,1)==-1  % ����߽�
    BCType=[BCType 'q'];
    Bifur(i,1:2)=2;
  else
    BCType=[BCType 'B'];  % �������嶯��ѭ�����ԣ�ֻ��B��û��C
    Bifur(i,1:2)=ConMat(i,1:2);
  end
  
  if ConMat(i,3)==-2  % ����߽�
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

%% %% ����Ѫ�ܲ��� %%%%
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
% ����1Dģ�͵�ճ����Ѫ�ܱ�ģ�Ͳ���
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

%% ͳ��R,C,L��ֵ
R=8*VesParam(7,:).*VesParam(1,:)./(pi.*(VesParam(2,:)/2).^4);
C=3*pi*(VesParam(2,:)/2).^3.*VesParam(1,:)./(2.*Eh');

%% ���ɱ߽����������ļ� %%%%
% ģ�ʹ������ɵ��ļ��ж�ȡ�߽粨��
inFileName=[NetTypeName '_IN'];
i_artinput=GenBoundInput_Systemic(inFileName,ModelParam);

fileName=GenInputFile_Systemic(NetTypeName,VesParam,BCTypeAll,Bifur,BCVal,ModelParam);
system(['mkdir C:\1DWin\' NetTypeName]);
system(['copy ' NetTypeName '*.* C:\1DWin\' NetTypeName '\']);

if RunFreqAnal==1
  AllNodes=union(From,To);
  JNum=length(AllNodes);
  df=0.5;    % Ƶ�ʲ���
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
    
    % ����ÿ��Ѫ�ܵ�ƽ��Ѫѹ����βѹ��
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
    % Norm�����壺��0HzΪ���գ�����Ƶ�ʵ����ߣ�˥����α仯
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
