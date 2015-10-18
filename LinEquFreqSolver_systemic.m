function X = LinEquFreqSolver_systemic(s,R,C,VesNum,JNum,AllNodes,From,To,OutR,OutC)
% 初始化AX=B中的A，B
cnt=0;
JMat=zeros(JNum+VesNum,JNum+VesNum);  % 左边矩阵
RHS=zeros(JNum+VesNum,1); % 右边矩阵，顺序为Q1->QVesNum,P1->PJNum
for i=1:VesNum
  FromNode=find(AllNodes==From(i));
  ToNode=find(AllNodes==To(i));
  
  % 处理RQ2(s)=P1(s)-P2(s)
  % 流出边界都作为阻力力边界处理
  outInd1=find(To(i)==From);     % 查询从该段血管的终点流出的血管
  inInd1=find(From(i)==From);    % 查询其他从该段血管起点出发的血管
  inInd2=find(From(i)==To);      % 查询流入该段血管起点的血管
  if length(outInd1)==0          % 流出边界
    RL=OutR(i)*1e10;
    JMat(i,i)=R(i)+RL;
    JMat(i,FromNode+VesNum)=-1;
  elseif length(inInd2)==0  % 流入边界
    JMat(i,i)=R(i);
    JMat(i,FromNode+VesNum)=-1;
    JMat(i,ToNode+VesNum)=1;
  else
    JMat(i,i)=R(i);
    JMat(i,FromNode+VesNum)=-1;
    JMat(i,ToNode+VesNum)=1;
  end
  
  % 处理CsP1(s)=Q1(s)-Q2(s)
  % 分析入口节点
  if length(inInd1)==2
    % 如果有2条血管具有该起点，说明该节点是一个分叉节点
    % 此时必有且仅有一条血管流入该点，即inInd2
    inInd1(inInd1==i)=[]; % 删除本身这条血管的编号
    Cap=C(i)+C(inInd1); % C永远取后级血管的值
    JMat(FromNode+VesNum,FromNode+VesNum)=Cap*s;
    JMat(FromNode+VesNum,i)=1;
    JMat(FromNode+VesNum,inInd1)=1;
    JMat(FromNode+VesNum,inInd2)=-1;
  elseif length(inInd2)==2
    % 如果有2条血管流入该起点，说明该节点是一个汇聚节点
    Cap=C(i);
    JMat(FromNode+VesNum,FromNode+VesNum)=Cap*s;
    JMat(FromNode+VesNum,i)=1;
    JMat(FromNode+VesNum,inInd2(1))=-1;
    JMat(FromNode+VesNum,inInd2(2))=-1;
  elseif length(inInd2)==0
    % 如果没有血管流入该起点，说明该段血管是一个流入边界
    cnt=cnt+1;
    Cap=C(i);
    JMat(FromNode+VesNum,FromNode+VesNum)=Cap*s;
    JMat(FromNode+VesNum,i)=1;
    % 考虑叠加定理，其余输入为开路(TODO:未实现)
    RHS(FromNode+VesNum)=93e-6;
  else
    % 以上情况均不符，则为接合血管
    Cap=C(i);
    JMat(FromNode+VesNum,FromNode+VesNum)=Cap*s;
    JMat(FromNode+VesNum,i)=1;
    JMat(FromNode+VesNum,inInd2)=-1;
  end
  
  % 分析出口节点
  outInd1=find(To(i)==From); % 查询其他流入该段血管终点的血管
  outInd2=find(To(i)==To);   % 查询流出该段血管终点的血管
  if length(outInd1)==2
    % 如果该段血管流向两根血管，说明该节点是一个分叉节点
    Cap=C(outInd1(1))+C(outInd1(2));
    JMat(ToNode+VesNum,ToNode+VesNum)=Cap*s;
    JMat(ToNode+VesNum,i)=-1;
    JMat(ToNode+VesNum,outInd1(1))=1;
    JMat(ToNode+VesNum,outInd1(2))=1;
  elseif length(outInd2)==2
    % 如果有2条血管到达该终点，说明该节点是一个汇聚节点
    outInd2(outInd2==i)=[];
    Cap=C(outInd1);
    JMat(ToNode+VesNum,ToNode+VesNum)=Cap*s;
    JMat(ToNode+VesNum,i)=-1;
    JMat(ToNode+VesNum,outInd2)=-1;
    JMat(ToNode+VesNum,outInd1)=1;
  elseif length(outInd1)==0
    % 如果该终点不是任何血管的起点，说明该节点是一个流出边界
    % TODO(panqing): 只有主输出是这样设置的，其他输出边界应同输入边界一样对待
    if 1
      RL=OutR(i)*1e10;
      JMat(ToNode+VesNum,ToNode+VesNum)=1;
      JMat(ToNode+VesNum,outInd2)=-RL;
    else
      JMat(ToNode+VesNum,ToNode+VesNum)=1;
      JMat(ToNode+VesNum,outInd2)=-1;
    end
  else
    % 以上情况均不符，则为接合血管
    Cap=C(outInd2);
    JMat(ToNode+VesNum,ToNode+VesNum)=Cap*s;
    JMat(ToNode+VesNum,i)=-1;
    JMat(ToNode+VesNum,outInd2)=1;
  end
end
% disp(cnt);
X=JMat\RHS;
end
