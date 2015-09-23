function [InvSeg,DiffSeg,AllBoundFlow] = BoundOptPulsatile(NetTypeID,SegType,VesNum,From,To,Diam,Len,Visc,WallTh,E,Boundary,SS_Press,SS_Flow,SS_DeltaP)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit
%%%% �߽�����
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% �����������±꣬������߽粻�����Ż�
switch NetTypeID
  case {Net_546_Meas_ID,Net_546_ID}
    MainBoundInd=[23 28];
  case Egg_818_ID
    %         MainInputInd
end

% Ƶ�ʷ�������
df=2;    % ����
f=0:df:10;

% PSO�㷨��������
nItNum=7;                                   % ��������
nPSO=10;                                    % ����Ⱥ�����Ӹ���
nParDim=length(BoundFlow);                  % ÿ�����ӵ�ά��. �˴�����Ϊlength(BoundFlow),��ʵ���Ż�ʱ��ȥ��������������
AllBoundFlow=zeros(nParDim,nPSO,nItNum);    % ����λ��xi(�����ӵ�ֵ)
v=zeros(nParDim,nPSO,nItNum);               % �����ٶ�vi
AllBoundFlowLimit=zeros(nParDim,2);         % �����У�ÿ��Ԫ�ص����ֵ��Сֵ����ֹ��������仯���󣬻�ı䷽��
InvSeg=zeros(VesNum,length(f),nPSO,nItNum); % ÿ�ε�����ÿ��Ƶ�ʷ���������Ⱥ��ÿ�����������£�����Ѫ�ܱ��
DiffSeg=1e100*ones(nPSO,nItNum);
Pg=zeros(nParDim,1);                        % ĳ������ȫ������
Pi=zeros(nParDim,nPSO);                     % ���������Ρ���ʷ����

% ��ʼ������Ⱥ
% ��һ�ε������������Ⱥ�У�ÿ�����ӵ�ÿ��Ԫ�س�ʼ��Ϊ�����
for i=1:nPSO
  for j=1:length(BoundFlow)
    ind=find(j==MainBoundInd);
    if ~isempty(ind)
      AllBoundFlow(MainBoundInd(ind),i,1)=BoundFlow(MainBoundInd(ind));
    else
      % ���ڷ������롢������߽磬���������ʼֵ
      % TODO: ����裿Ŀǰ��Ϊԭ�߽�ֵ*�����*10
      AllBoundFlow(j,i,1)=BoundFlow(j)*rand;
    end
  end
end
% ����������Ԫ�ص������ޣ��߽��Ѫ�����򱣳ֲ���
AllBoundFlowLimit(:,1)=mean(AllBoundFlow(:,:,1),2)*1e-3;
AllBoundFlowLimit(:,2)=mean(AllBoundFlow(:,:,1),2)*5;
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

% ��λ����
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s
WallTh=WallTh*1e-6; %m

% ȥ���٣���֤���Է������Ч��
% m_ratio=1e6;
% kg_ratio=1e2;
m_ratio=1;
kg_ratio=1;
Len=Len*m_ratio;
Diam=Diam*m_ratio;
Visc=Visc*kg_ratio;
WallTh=WallTh*m_ratio;
E=E*kg_ratio;

% ׼��������Է�����
DeltaPX=zeros(VesNum,1);PressXX=DeltaPX;RatioPX=DeltaPX;RatioP2=DeltaPX;
% JMat=zeros(JNum+VesNum,JNum+VesNum);  % ��߾���
% RHS=zeros(JNum+VesNum,1); % �ұ߾���˳��ΪQ1->QVesNum,P1->PJNum
AllNodes=union(From,To);
JNum=length(AllNodes);

% �ٶȡ�λ�ø��¹�ʽ����
c1=2; % ��֪
c2=2; % ���
r=1;  % Լ������
w=0.729;  % ����Ȩ��

for m=1:nItNum
  fprintf(1,'ItNum=%d\n',m);
%   c1=c1*0.8;
%   c2=c2*1.2;
  for i=1:nPSO
    fprintf(1,'  nPSO=%d\n',i);
    % ��AllBoundFlowΪ�߽磬������̬Ѫ������ѧ���ԣ��Ӷ��õ���������߽��
    % SS_Press,SS_Flow��SS_DeltaP
    Boundary(:,3)=AllBoundFlow(:,i,m);
    [SS_Press,SS_Flow,SS_Vel,SS_DeltaP,SS_Visc,SS_Hd]=...
      LinEqu(NetTypeID,DampInit,VesNum,SegType,From,To,Diam,Len,Visc,Boundary);
    % ����������˳Ӧ�Բ���
    R=(8.*Visc.*Len)./(pi*(Diam./2).^4);
    C=3*pi*(Diam./2).^3.*Len./(2.*E.*WallTh);
    AllPressX=zeros(length(f),VesNum);AllFlowX=AllPressX;AllDeltaPX=AllPressX;AllRatioPX=AllPressX;
    for k=1:length(f)
      fprintf(1,'    Freq=%.3fHz\n',f(k));
      s=1i*2*pi*f(k);
      % ע��LinEquFreqSolver������Ȼ������SS_Press,SS_Flow��SS_DeltaP����ʵ����ֻ�õ����⼸�������еı߽�Ѫ��ֵ
      X=LinEquFreqSolver(s,R,C,VesNum,JNum,AllNodes,From,To,SS_Press,SS_Flow,SS_DeltaP,m_ratio,kg_ratio);
      absX=abs(X);  % ȡģ
      FlowX=absX(1:VesNum);
      PressX=absX(VesNum+1:end);
      
      % ����ÿ��Ѫ�ܵ�ƽ��Ѫѹ����βѹ��
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
      % �涨������Ϊ1������Ϊ0
      InvInd=find(AllDeltaPX(k,:)<0); % �ҵ���k��Ƶ�ʵķ������У�Ѫ������Ϊ����Ѫ�ܶα��
      InvSeg(InvInd,k,i,m)=1;         % ��¼��k��Ƶ�ʵķ�������������Ѫ��
    end
    % ����Ƶ�������Ҫ�ѽ�������һ��ֵ��DiffSeg(nPSO,nItNum)
%     DiffSeg(i,m)=mean(std(InvSeg(:,:,i,m),0,2));  % ����ÿ��Ƶ���£�Ѫ������ı�����
    DiffSeg(i,m)=sum(sum(InvSeg(:,:,i,m),2));
  end
  % �Ż�Ŀ��Ϊ����Ѫ�ܶ���
  [Value1,Ind1]=min(DiffSeg);          % IndΪȫ����������������ʷ����ֵPg�����
  [Value2,Ind2]=min(Value1);
  Pg=AllBoundFlow(:,Ind1(Ind2),Ind2);    % PgΪ����(m)�����е�ȫ������ֵ���Ƿ�Ӧ��ȡ��ʷ�����е�ȫ������ֵ����
  % Ѱ����ʷ��ÿ�����ӵ�����ֵ
  for k=1:nPSO
    [Value,Ind]=min(DiffSeg(k,:));   % IndΪ��k�����ӵ���ʷ����ֵ�����
    Pi(:,k)=AllBoundFlow(:,k,Ind);  % Pi(:,k)Ϊ��k�����ӵ���ʷ����ֵ
    v(:,k,m+1)=w.*v(:,k,m)+c1*rand*(Pi(:,k)-AllBoundFlow(:,k,m))+c2*rand*(Pg-AllBoundFlow(:,k,m));
    AllBoundFlow(:,k,m+1)=AllBoundFlow(:,k,m)+r*v(:,k,m+1);
    % ���Ʊ߽������䶯��Χ
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
