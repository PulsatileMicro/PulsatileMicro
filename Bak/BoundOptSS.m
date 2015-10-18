function [InvSeg,AllBoundFlow] = BoundOptSS(DampFactor,NetTypeID,VesNum,SegType,From,To,Diam,Len,Visc,Boundary)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN

% �߽�����
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% PSO�㷨��������
nItNum=20;                    % ��������
nPSO=25;                      % ����Ⱥ�����Ӹ���
nParDim=length(BoundFlow);    % ÿ�����ӵ�ά��. �˴�����Ϊlength(BoundFlow),��ʵ���Ż�ʱ��ȥ��������������
AllBoundFlow=zeros(nParDim,nPSO,nItNum);  % ����λ��xi
v=zeros(nParDim,nPSO,nItNum);             % �����ٶ�vi
AllBoundFlowLimit=zeros(nParDim,2);       % �����У�ÿ��Ԫ�ص����ֵ��Сֵ
InvSeg=VesNum*ones(nPSO,nItNum);  % ��Ӧ�Ƚ����������Ѫ�ܶ�������ֵΪVesNum��
InvSegVel=1e100*ones(nPSO,nItNum);% ������Ӧ�Ƚ����������Ѫ�ܵ����ٺ�
Pg=zeros(nParDim,1);              % ĳ������ȫ������
Pi=zeros(nParDim,nPSO);           % ��������ʷ����
AllMeanP=zeros(VesNum,nItNum);
AllMeanFlow=zeros(VesNum,nItNum);

% �߽�Ѫ�ܷ������� %%%%
% ���ݱ߽������ķ�������
% Ѫ��Ϊ����->��߽�->�߽�ڵ���From����
% Ѫ��Ϊ����->���߽�->�߽�ڵ���To����
[From,To]=OriNodeModify(Boundary,From,To);
Qref=0.001; % TODO: Qrefû�ã����������Ҫ�ع�
% Boundary(:,4)=0.45;
BHd=BoundaryInput(Boundary,From,To,Qref);

% ��λ����
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s

% �����������±�(�ڱ߽������е��±�)��������߽粻�����Ż�
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
% ��ʼ������Ⱥ
% ÿ�������е�ÿ��Ԫ�س�ʼ��Ϊ�����
for i=1:nPSO
  for j=1:length(BoundFlow)
    ind=find(j==MainBoundInd);
    if ~isempty(ind)
      AllBoundFlow(MainBoundInd(ind),i,1)=BoundFlow(MainBoundInd(ind));
    else
      % ���ڷ������롢������߽磬���������ʼֵ
      % TODO: ����裿
      AllBoundFlow(j,i,1)=BoundFlow(j)*rand*1e2;
    end
  end
end
% ����������Ԫ�ص������ޣ��߽��Ѫ�����򱣳ֲ���
AllBoundFlowLimit(:,1)=mean(AllBoundFlow(:,:,1),2)*1e-3;
AllBoundFlowLimit(:,2)=mean(AllBoundFlow(:,:,1),2)*1e3;
AllBoundFlowLimit(MainBoundInd,1)=AllBoundFlow(MainBoundInd,1,1);
AllBoundFlowLimit(MainBoundInd,2)=AllBoundFlow(MainBoundInd,1,1);
% ����߽�Ϊ������Ԫ�������޶Ե���ȷ���Ƚϴ�Сʱ�߼���ȷ
for i=1:length(BoundFlow)
  if BoundFlow(i)<0
    temp=AllBoundFlowLimit(i,1);
    AllBoundFlowLimit(i,1)=AllBoundFlowLimit(i,2);
    AllBoundFlowLimit(i,2)=temp;
  end
end
% �ٶȡ�λ�ø��¹�ʽ����
c1=4; % ��֪
c2=2; % ���
r=1;  % Լ������
w=0.729;  % ����Ȩ��
CTime=5;
for m=1:nItNum  % ����nItNum��
  m
  c1=c1*0.8;
  c2=c2*1.2;
  for k=1:nPSO  % ����ÿһ������
    T=0;
    Visc=2*ones(VesNum,1)*1e-3;       % ��ʼճ�Ͷ�
    Porder=1:VesNum;Norder=VesNum:1;  % ��ʼHd����˳��
    errorFlag=0;
    while T<CTime   % ѭ��1��Visc����
      T=T+1
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
      
      if NetTypeID~=Net_546_ID && NetTypeID~=Net_122_ID
        %%ճ�Ͷȷ���
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
      
      % ��¼ÿ�ε����ķ�����
      DebugVisc(:,T)=Visc;
      DebugHd(:,T)=Hd;
      DebugPressure(:,T)=MeanP;
      DebugFlow(:,T)=MeanFlowNew;
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
end