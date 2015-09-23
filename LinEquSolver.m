function [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam,Visc,Len,From,To,~)
J=(pi.*(Diam./2).^4)./(8.*Visc.*Len);   % 546段血管的conductance--J, 依据数据文件顺序

AllNodes=union(From,To);
BifurOutFlag=0; % 主输出为分叉的标志
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
% 初始化AX=B的A和B分别为JMat和NodeFlow
JMat=zeros(JNum,JNum);  % 左边矩阵
NodeFlow=zeros(JNum,1); % 右边矩阵
for i=1:length(AllNodes)
  BoundNodeIndex=find(AllNodes(i)==BoundNode);    % 寻找节点是否属于边界节点。若是，BoundNodeIndex为边界矩阵中的序号
  OutSegIndex=find(From==AllNodes(i));         % 当前节点为一段血管的流出节点时，血管段的序号（From向量的序号==血管段的序号）
  InSegIndex=find(To==AllNodes(i));          % 当前节点为一段血管的流入节点时，血管段的序号
  FromNode=[];ToNode=[];BifurOutFlag=0;
  for j=1:length(OutSegIndex)
    Ind=find(To(OutSegIndex(j))==AllNodes);
    if ~isempty(Ind)
      ToNode(j)=Ind;     % ToNode记录该血管段流向的另一端节点的序号
    end
  end
  for j=1:length(InSegIndex)
    Ind=find(From(InSegIndex(j))==AllNodes);
    if ~isempty(Ind)
      FromNode(j)=Ind;  % FromNode记录流向该节点的血管的另一端节点的序号
    end
  end
  
  if ~isempty(BoundNodeIndex)                % 边界节点
    if BoundType(BoundNodeIndex)~=0          % 流量节点
      % 流入流量节点
      if ~isempty(OutSegIndex)
        JMat(i,ToNode)=-J(OutSegIndex);
        JMat(i,i)=J(OutSegIndex);
        NodeFlow(i)=BoundFlow(BoundNodeIndex)/1e12/60;
      % 流出流量节点
      elseif ~isempty(InSegIndex)
        JMat(i,FromNode)=-J(InSegIndex);
        JMat(i,i)=J(InSegIndex);
        NodeFlow(i)=BoundFlow(BoundNodeIndex)/1e12/60;
      else
        error('error!');
      end
    end
  else                              % 非边界节点
    % 找到压力输出边界前的一点
    POutInd=find(To==POutNode);
    PPreOutNode=From(POutInd);
    if AllNodes(i)==PPreOutNode
      if length(InSegIndex)==1  % 一进两出的情况，CAM网络中出现
        BifurOutFlag=1;   % 主输出为分叉的标志
        NodeFlow(i)=-POut*J(POutInd)*133;
      else
        BifurOutFlag=0;
        NodeFlow(i)=-POut*J(OutSegIndex)*133;
      end
    end
    if length(OutSegIndex)==2   % 一进二出的情况，入节点1段，出节点2段
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
        % 压力边界情况
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
    elseif length(InSegIndex)==2 % 二进一出的情况，入节点2段，出节点1段
      JMat(i,i)=-(J(InSegIndex(1))+J(InSegIndex(2))+J(OutSegIndex));
      JMat(i,ToNode)=J(OutSegIndex);
      % 同上，二进一出的二进可能来自同一点，此时需将两个J相加，避免J被覆盖掉
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
    elseif length(InSegIndex)==3  %三进零出
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
    elseif length(OutSegIndex)==3  %零进三出
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
    else  % 一进一出的情况，入节点1段，出节点1段
      JMat(i,i)=-(J(InSegIndex)+J(OutSegIndex));
      JMat(i,FromNode)=J(InSegIndex);
      JMat(i,ToNode)=J(OutSegIndex);
    end
  end
end

P=JMat\NodeFlow;

% 依据P是否为实数判断结果是否有效
if isreal(P)
  Eju=0;
else
  Eju=1;
end

MeanP=zeros(VesNum,1);
MeanFlow=zeros(VesNum,1);
DeltaP=zeros(VesNum,1);

% 计算每段血管的平均血压和首尾压差
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

% 计算每段血管的流量
for j=1:VesNum
  MeanFlow(j)=1e12*60*J(j)*DeltaP(j);  %nl/min
end

DeltaP=DeltaP/133;    %单位转换mmHg
MeanP=MeanP/133;      %单位转换mmHg