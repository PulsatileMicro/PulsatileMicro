function [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam,Visc,Len,From,To,~)
J=(pi.*(Diam./2).^4)./(8.*Visc.*Len);   % 546��Ѫ�ܵ�conductance--J, ���������ļ�˳��

AllNodes=union(From,To);
BifurOutFlag=0; % �����Ϊ�ֲ�ı�־
for i=1:length(BoundNode)
  if BoundType(i)==0
    POutInd=find(AllNodes==BoundNode(i));
    if POutInd
      POutNode=AllNodes(POutInd);
      POut=BoundFlow(i);
      AllNodes(POutInd)=[];
    end
  end
end
JNum=length(AllNodes);
% ��ʼ��AX=B��A��B�ֱ�ΪJMat��NodeFlow
JMat=zeros(JNum,JNum);  % ��߾���
NodeFlow=zeros(JNum,1); % �ұ߾���
for i=1:length(AllNodes)
  BoundNodeIndex=find(AllNodes(i)==BoundNode);    % Ѱ�ҽڵ��Ƿ����ڱ߽�ڵ㡣���ǣ�BoundNodeIndexΪ�߽�����е����
  OutSegIndex=find(From==AllNodes(i));         % ��ǰ�ڵ�Ϊһ��Ѫ�ܵ������ڵ�ʱ��Ѫ�ܶε���ţ�From���������==Ѫ�ܶε���ţ�
  InSegIndex=find(To==AllNodes(i));          % ��ǰ�ڵ�Ϊһ��Ѫ�ܵ�����ڵ�ʱ��Ѫ�ܶε����
  FromNode=[];ToNode=[];BifurOutFlag=0;
  for j=1:length(OutSegIndex)
    Ind=find(To(OutSegIndex(j))==AllNodes);
    if ~isempty(Ind)
      ToNode(j)=Ind;     % ToNode��¼��Ѫ�ܶ��������һ�˽ڵ�����
    end
  end
  for j=1:length(InSegIndex)
    Ind=find(From(InSegIndex(j))==AllNodes);
    if ~isempty(Ind)
      FromNode(j)=Ind;  % FromNode��¼����ýڵ��Ѫ�ܵ���һ�˽ڵ�����
    end
  end
  
  if ~isempty(BoundNodeIndex)                % �߽�ڵ�
    if BoundType(BoundNodeIndex)~=0          % �����ڵ�
      % ���������ڵ�
      if ~isempty(OutSegIndex)
        JMat(i,ToNode)=-J(OutSegIndex);
        JMat(i,i)=J(OutSegIndex);
        NodeFlow(i)=BoundFlow(BoundNodeIndex)/1e12/60;
      % ���������ڵ�
      elseif ~isempty(InSegIndex)
        JMat(i,FromNode)=-J(InSegIndex);
        JMat(i,i)=J(InSegIndex);
        NodeFlow(i)=BoundFlow(BoundNodeIndex)/1e12/60;
      else
        error('error!');
      end
    end
  else                              % �Ǳ߽�ڵ�
    % �ҵ�ѹ������߽�ǰ��һ��
    POutInd=find(To==POutNode);
    PPreOutNode=From(POutInd);
    if AllNodes(i)==PPreOutNode
      if length(InSegIndex)==1  % һ�������������CAM�����г���
        BifurOutFlag=1;   % �����Ϊ�ֲ�ı�־
        NodeFlow(i)=-POut*J(POutInd)*133;
      else
        BifurOutFlag=0;
        NodeFlow(i)=-POut*J(OutSegIndex)*133;
      end
    end
    if length(OutSegIndex)==2   % һ���������������ڵ�1�Σ����ڵ�2��
      JMat(i,i)=-(J(InSegIndex)+J(OutSegIndex(1))+J(OutSegIndex(2)));
      JMat(i,FromNode)=J(InSegIndex);
      if(ToNode(1)==ToNode(2))
        JMat(i,ToNode(1))=J(OutSegIndex(1))+J(OutSegIndex(2));
      elseif (ToNode(1)==FromNode)
        JMat(i,FromNode)=J(InSegIndex)+J(OutSegIndex(1));
        JMat(i,ToNode(2))=J(OutSegIndex(2));
      elseif (ToNode(2)==FromNode)
        JMat(i,FromNode)=J(InSegIndex)+J(OutSegIndex(2));
        JMat(i,ToNode(1))=J(OutSegIndex(1));
      else
        % ѹ���߽����
        if BifurOutFlag
          DelInd=(OutSegIndex==POutInd);
          OutSegIndex(DelInd)=[];
          if ToNode(1)>0
            JMat(i,ToNode(1))=J(OutSegIndex);
          else
            JMat(i,ToNode(2))=J(OutSegIndex);
          end
        else
          JMat(i,ToNode(1))=J(OutSegIndex(1));
          JMat(i,ToNode(2))=J(OutSegIndex(2));
        end
      end
    elseif length(InSegIndex)==2 % ����һ�����������ڵ�2�Σ����ڵ�1��
      JMat(i,i)=-(J(InSegIndex(1))+J(InSegIndex(2))+J(OutSegIndex));
      JMat(i,ToNode)=J(OutSegIndex);
      % ͬ�ϣ�����һ���Ķ�����������ͬһ�㣬��ʱ�轫����J��ӣ�����J�����ǵ�
      if(FromNode(1)==FromNode(2))
        JMat(i,FromNode(1))=J(InSegIndex(1))+J(InSegIndex(2));
      elseif (FromNode(1)==ToNode)
        JMat(i,ToNode)=J(OutSegIndex)+J(InSegIndex(1));
        JMat(i,FromNode(2))=J(InSegIndex(2));
      elseif (FromNode(2)==ToNode)
        JMat(i,ToNode)=J(OutSegIndex)+J(InSegIndex(2));
        JMat(i,FromNode(1))=J(InSegIndex(1));
      else
        JMat(i,FromNode(1))=J(InSegIndex(1));
        JMat(i,FromNode(2))=J(InSegIndex(2));
      end
    elseif length(InSegIndex)==3  %�������
      JMat(i,i)=-(J(InSegIndex(1))+J(InSegIndex(2))+J(InSegIndex(3)));
      if (FromNode(1)==FromNode(2))
        JMat(i,FromNode(1))=J(InSegIndex(1))+J(InSegIndex(2));
        JMat(i,FromNode(3))=J(InSegIndex(3));
      elseif (FromNode(2)==FromNode(3))
        JMat(i,FromNode(2))=J(InSegIndex(2))+J(InSegIndex(3));
        JMat(i,FromNode(1))=J(InSegIndex(1));
      elseif (FromNode(1)==FromNode(3))
        JMat(i,FromNode(1))=J(InSegIndex(1))+J(InSegIndex(3));
        JMat(i,FromNode(2))=J(InSegIndex(2));
      else
        JMat(i,FromNode(1))=J(InSegIndex(1));
        JMat(i,FromNode(2))=J(InSegIndex(2));
        JMat(i,FromNode(3))=J(InSegIndex(3));
      end
    elseif length(OutSegIndex)==3  %�������
      JMat(i,i)=-(J(OutSegIndex(1))+J(OutSegIndex(2))+J(OutSegIndex(3)));
      if (ToNode(1)==ToNode(2))
        JMat(i,ToNode(1))=J(OutSegIndex(1))+J(OutSegIndex(2));
        JMat(i,ToNode(3))=J(OutSegIndex(3));
      elseif (ToNode(2)==ToNode(3))
        JMat(i,ToNode(2))=J(OutSegIndex(2))+J(OutSegIndex(3));
        JMat(i,ToNode(1))=J(OutSegIndex(1));
      elseif (ToNode(1)==ToNode(3))
        JMat(i,ToNode(1))=J(OutSegIndex(1))+J(OutSegIndex(3));
        JMat(i,ToNode(2))=J(OutSegIndex(2));
      else
        JMat(i,ToNode(1))=J(OutSegIndex(1));
        JMat(i,ToNode(2))=J(OutSegIndex(2));
        JMat(i,ToNode(3))=J(OutSegIndex(3));
      end
    else  % һ��һ�����������ڵ�1�Σ����ڵ�1��
      JMat(i,i)=-(J(InSegIndex)+J(OutSegIndex));
      JMat(i,FromNode)=J(InSegIndex);
      JMat(i,ToNode)=J(OutSegIndex);
    end
  end
end

P=JMat\NodeFlow;

% ����P�Ƿ�Ϊʵ���жϽ���Ƿ���Ч
if isreal(P)
  Eju=0;
else
  Eju=1;
end

MeanP=zeros(VesNum,1);
MeanFlow=zeros(VesNum,1);
DeltaP=zeros(VesNum,1);

% ����ÿ��Ѫ�ܵ�ƽ��Ѫѹ����βѹ��
for j=1:VesNum
  InIndex=find(From(j)==AllNodes);
  OutIndex=find(To(j)==AllNodes);
  if ~isempty(OutIndex)
    DeltaP(j)=P(InIndex)-P(OutIndex);
    MeanP(j)=P(InIndex)-DeltaP(j)/2;
  else
    OutInd=find(To(j)==BoundNode);
    DeltaP(j)=P(InIndex)-BoundFlow(OutInd)*133;
    MeanP(j)=P(InIndex)-DeltaP(j)/2;
  end
end

% ����ÿ��Ѫ�ܵ�����
for j=1:VesNum
  MeanFlow(j)=1e12*60*J(j)*DeltaP(j);  %nl/min
end

DeltaP=DeltaP/133;    %��λת��mmHg
MeanP=MeanP/133;      %��λת��mmHg