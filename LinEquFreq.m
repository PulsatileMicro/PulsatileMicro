function [AllFlowX,AllPressX,NormAllFlowX,NormAllPressX,CutFreqFlow,CutFreqPress,AllDeltaPX,AllRatioPX]=...
  LinEquFreq(VesNum,From,To,Diam,Len,Visc,WallTh,E,Boundary,SS_Press,SS_Flow,SS_DeltaP,DampFactor,ViscRatio,InputType,f,df)
%% Ѫ�������Ƶ�����
%%%% �߽�����
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% ��λ����
Len=Len*1e-6;       %m
Diam=Diam*1e-6;     %m
Visc=Visc*1e-3;     %Pa.s
WallTh=WallTh*1e-6; %m

% ȥ���٣���֤���Է������Ч��(TODO: δʵ��)
len_ratio=1;
mass_ratio=1;
Len=Len*len_ratio;
Diam=Diam*len_ratio;
Visc=Visc*mass_ratio;
WallTh=WallTh*len_ratio;
E=E*mass_ratio;

% ����������˳Ӧ��
R=(8.*Visc.*Len)./(pi*(Diam./2).^4);
C=3*pi*(Diam./2).^3.*Len./(2.*E.*WallTh);

DeltaPX=zeros(VesNum,1);PressXX=DeltaPX;RatioPX=DeltaPX;RatioP2=DeltaPX;
% ׼��������Է�����
% AllNodes=union(From,To);
% JNum=length(AllNodes);
% JMat=zeros(JNum+VesNum,JNum+VesNum);  % ��߾���
% RHS=zeros(JNum+VesNum,1); % �ұ߾���˳��ΪQ1->QVesNum,P1->PJNum
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
  
  % ����ÿ��Ѫ�ܵ�ƽ��Ѫѹ����βѹ��
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

% ����ÿ��Ѫ�ܵĽ�ֹƵ��
CutFreqFlow=zeros(VesNum,1);CutFreqPress=zeros(VesNum,1);
% ÿ��Ƶ�ʷ����ϵ�ֵ��������Ѫ������һ��
% TODO: ���Ѫ�ܲ�һ���Ǳ��1
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

