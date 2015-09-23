%% �Ż�Ѫ������߽�����
function [AllBoundFlow,InvSeg] = BOUND_OPT_FUNC(NetTypeID,Boundary,DatMatrix)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global NOT_OPT_PARA OPT_PARA OPT_BOUND_FLOW OPT_BOUND_META
OptType=OPT_BOUND_FLOW;
%% 1. ���ݲ�ͬ���͵�Ѫ������ȷ�������롢������߽�
% ��ȷ���ı��Ϊ����ţ�����Ѫ�ܶ�����
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

%% 2. ���ݽ��
% Ѫ�ܶ�����
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

% �߽�����
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
% BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);
BoundMainInOut=zeros(length(BoundNode),1);

for i=1:length(BoundNode)
  % ���ұ߽�ڵ�Ϊ��ڽڵ��Ѫ�ܵı�ţ�Ϊ���ҵ�������
  SegInd=find(From==BoundNode(i));
  if BoundType(i)==0
    % ѹ���߽�
    BoundFlow(i)=Boundary(i,3);
    BoundMainInOut(i)=1;
  elseif ~isempty(SegInd)
    % �����ýڵ�Ϊ����ڵ�
    if ~isempty(find(MainIn==SegInd,1))
      % �鵽������
      BoundFlow(i)=Boundary(i,3);
      BoundMainInOut(i)=1;
    else
      % һ������
%       BoundFlow(i)=(0.4*Diam(SegInd)-1.9)*0.25*pi.*Diam(SegInd).^2/1e3*60;
      BoundFlow(i)=30;
    end
  else
    SegInd=find(To==BoundNode(i));
    %  ����ڵ�
%     BoundFlow(i)=-(0.4*Diam(SegInd)-1.9)*0.25*pi.*Diam(SegInd).^2/1e3*60;
    BoundFlow(i)=-30;  % nL/min
  end
end

% ��λ����
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s

% �߽�Ѫ�ܷ������� %%%%
% ���ݱ߽������ķ�������
% Ѫ��Ϊ����->��߽�->�߽�ڵ���From����
% Ѫ��Ϊ����->���߽�->�߽�ڵ���To����
[From,To]=OriNodeModify(Boundary,From,To);

Qref=0.001; % TODO: Qrefû�ã����������Ҫ�ع�
% Boundary(:,4)=0.45;
BHd=BoundaryInput(Boundary,From,To,Qref);

%% ճ�Ͷ�ѭ��
% ѭ������
Loop1_Cnt=0;
Loop1_Visc_Limit=20;

DebugVisc=zeros(VesNum,Loop1_Visc_Limit);
DebugHd=zeros(VesNum,Loop1_Visc_Limit);
DebugPressure=zeros(VesNum,Loop1_Visc_Limit);
DebugFlow=zeros(VesNum,Loop1_Visc_Limit);
Porder=1:VesNum;Norder=VesNum:1;  % ��ʼHd����˳��

%% ����Ⱥ�㷨����
nItNum=20;                    % ��������
nPSO=25;                      % ����Ⱥ�����Ӹ���
nParDim=length(BoundFlow);    % ÿ�����ӵ�ά��
AllBoundFlow=zeros(nParDim,nPSO,nItNum);  % ��������xi
AllBoundFlowLimit=zeros(nParDim,2);       % ���ӵ����ֵ��Сֵ
v=zeros(nParDim,nPSO,nItNum);             % �����ٶ�vi

for i=1:nPSO
  % ��ʼ������Ⱥ
  for j=1:length(BoundFlow)
    if BoundMainInOut(j)
      AllBoundFlow(j,i,1)=BoundFlow(j);
    else
      AllBoundFlow(j,i,1)=BoundFlow(j)*rand;
    end
  end
end
% ����������Ԫ�ص�������
% �߽��Ѫ�����򱣳ֲ���
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
c1=2; % ��֪
c2=2; % ���
r=1;  % Լ������
w=0.729;  % ����Ȩ��
InvSeg=VesNum*ones(nPSO,nItNum);  % ��Ӧ�Ƚ����������Ѫ�ܶ�������ֵΪVesNum��
InvSegVel=1e100*ones(nPSO,nItNum);% ������Ӧ�Ƚ����������Ѫ�ܵ����ٺ�
Pg=zeros(nParDim,1);              % ĳ������ȫ������
Pi=zeros(nParDim,nPSO);           % ��������ʷ����
AllMeanP=zeros(VesNum,nItNum);
AllMeanFlow=zeros(VesNum,nItNum);

for m=1:nItNum  % ����nItNum��
%   c1=c1*0.8;
%   c2=c2*1.2;
  for k=1:nPSO  % ����ÿһ������
    fprintf('Iter=%d, PSOid=%d ViscLoopCnt=%d\n',m,k,Loop1_Cnt);
    Loop1_Cnt=0;
    Loop1_Visc_Num=Loop1_Visc_Limit;
    Visc=2*ones(VesNum,1)*1e-3;
    Porder=1:VesNum;Norder=VesNum:1;  % ��ʼHd����˳��
    errorFlag=0;
    while Loop1_Cnt<Loop1_Visc_Num   % ѭ��1��Visc����
      Loop1_Cnt=Loop1_Cnt+1;
      % ���Է������ģ��
      [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,AllBoundFlow(:,k,m),VesNum,Diam,Visc,Len,From,To);
      % ��������ģ��
      [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju);
      %     FromNew=From;ToNew=To;MeanFlowNew=MeanFlow;DeltaPNew=DeltaP;
      % ����Hd����˳��
      [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,2,Eju);
      if Eju==1
        errorFlag=1;
        break;
      end
      % ����Hd
      [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam,MeanFlowNew,Eju);
      % ճ�Ͷȷ���
      umDiam=Diam.*1e6;   % um
      for i=1:VesNum
        Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
      end
      DebugVisc(:,Loop1_Cnt)=Visc;
      % ѭ������ڶ��κ�(Loop1_Cnt>1)����ʼУ��Visc
      % ����SOR����У�����������ճ�Ͷ��񵴵����
      ViscModType=1;  % 0����������1����Bifurcation������2����All segments����
      Alpha=0.5;
      DebugVisc=modifyViscosity(Porder,FromNew,ToNew,DebugVisc,Loop1_Cnt,ViscModType,Alpha,Eju);
      Visc=DebugVisc(:,Loop1_Cnt);
      
      %%%% �ж�ѭ������ģ�� %%%%
      % LoopOutType:0-������1-������2-�ﵽ����δ������3-visc����4-���Է��̴���
      [ViscMAE,Loop1OutType,Loop1_Visc_Num]=ViscLoopTerminator(DebugVisc*1e3,Loop1_Visc_Num,Loop1_Cnt,Eju);
    end
    if errorFlag~=1
      InvSeg(k,m)=length(find(MeanFlow<0));   % ÿһ�����Ӷ�Ӧ����Ӧ�Ⱥ����Ľ��
      MeanVel=MeanFlow./(0.25*pi*Diam.^2);
      InvSegVel(k,m)=-sum(MeanVel(MeanVel<0));
    end
  end
  % �Ż�Ŀ��Ϊ����Ѫ�ܶ���
  [Value1,Ind1]=min(InvSeg);          % IndΪȫ����������������ʷ����ֵPg�����
  [Value2,Ind2]=min(Value1);
  Pg=AllBoundFlow(:,Ind1(Ind2),Ind2);         % PgΪ����(m)�����е�ȫ������ֵ���Ƿ�Ӧ��ȡ��ʷ�����е�ȫ������ֵ����
  % Ѱ����ʷ��ÿ�����ӵ�����ֵ
  for k=1:nPSO
    [Value,Ind]=min(InvSeg(k,:));   % IndΪ��k�����ӵ���ʷ����ֵ�����
    Pi(:,k)=AllBoundFlow(:,k,Ind);  % Pi(:,k)Ϊ��k�����ӵ���ʷ����ֵ
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
  % ��¼���ε����е���СMeanP��MeanFlow
  AllMeanP(:,m)=MeanP;
  AllMeanFlow(:,m)=MeanFlow;
end

Vel=MeanFlowNew/1e12/60./(pi*0.25*Diam.^2)*1e3;