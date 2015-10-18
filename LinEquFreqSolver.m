function X=LinEquFreqSolver(s,R,C,VesNum,JNum,AllNodes,From,To,MeanP,Flow,DeltaP,m_ratio,kg_ratio,DampFactor,ViscRatio,InputType,Boundary)
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN E5 E02 Visc5 Visc02 DampE DampVisc
% ��ʼ��AX=B�е�A��B
cnt=0;
JMat=zeros(JNum+VesNum,JNum+VesNum);  % ��߾���
RHS=zeros(JNum+VesNum,1); % �ұ߾���˳��ΪQ1->QVesNum,P1->PJNum
for i=1:VesNum
  FromNode=find(AllNodes==From(i));
  ToNode=find(AllNodes==To(i));
  
  % ����RQ2(s)=P1(s)-P2(s)
  % �����߽綼��Ϊ�������߽紦��
  outInd1=find(To(i)==From);     % ��ѯ�Ӹö�Ѫ�ܵ��յ�������Ѫ��
  inInd1=find(From(i)==From);    % ��ѯ�����Ӹö�Ѫ����������Ѫ��
  inInd2=find(From(i)==To);      % ��ѯ����ö�Ѫ������Ѫ��
  if length(outInd1)==0          % �����߽�
    RL=MeanP(i)*133/kg_ratio/(Flow(i)/60/1e12*m_ratio^2);
    switch DampFactor
      case DampVisc
        RL=RL*ViscRatio;
    end
    JMat(i,i)=R(i)+RL;
    JMat(i,FromNode+VesNum)=-1;
  elseif length(inInd2)==0  % ����߽�
    if InputType==0   % ��������߽�
      JMat(i,i)=R(i);
      JMat(i,FromNode+VesNum)=-1;
      JMat(i,ToNode+VesNum)=1;
    else  % ѹ������߽�
      JMat(i,i)=R(i);
      % JMat(i,FromNode+VesNum)=-1;
      JMat(i,ToNode+VesNum)=1;
      RHS(i)=MeanP(i)*133/kg_ratio;
    end
  else
    JMat(i,i)=R(i);
    JMat(i,FromNode+VesNum)=-1;
    JMat(i,ToNode+VesNum)=1;
  end
  
  % ����CsP1(s)=Q1(s)-Q2(s)
  % ������ڽڵ�
  if length(inInd1)==2
    % �����2��Ѫ�ܾ��и���㣬˵���ýڵ���һ���ֲ�ڵ�
    % ��ʱ�����ҽ���һ��Ѫ������õ㣬��inInd2
    inInd1(inInd1==i)=[]; % ɾ����������Ѫ�ܵı��
    Cap=C(i)+C(inInd1); % C��Զȡ��Ѫ�ܵ�ֵ
    JMat(FromNode+VesNum,FromNode+VesNum)=Cap*s;
    JMat(FromNode+VesNum,i)=1;
    JMat(FromNode+VesNum,inInd1)=1;
    JMat(FromNode+VesNum,inInd2)=-1;
  elseif length(inInd2)==2
    % �����2��Ѫ���������㣬˵���ýڵ���һ����۽ڵ�
    Cap=C(i);
    JMat(FromNode+VesNum,FromNode+VesNum)=Cap*s;
    JMat(FromNode+VesNum,i)=1;
    JMat(FromNode+VesNum,inInd2(1))=-1;
    JMat(FromNode+VesNum,inInd2(2))=-1;
  elseif length(inInd2)==0
    % ���û��Ѫ���������㣬˵���ö�Ѫ����һ������߽�
    cnt=cnt+1;
    Cap=C(i);
    if InputType==0 % ��������߽�
      JMat(FromNode+VesNum,FromNode+VesNum)=Cap*s;
      JMat(FromNode+VesNum,i)=1;
      % ���ǵ��Ӷ�����������Ϊ��·(TODO:δʵ��)
      RHS(FromNode+VesNum)=Flow(i)/60/1e12*m_ratio^2;
    else  % ѹ������߽�
      JMat(FromNode+VesNum,FromNode+VesNum)=1;
      RHS(FromNode+VesNum)=MeanP(i)*133/kg_ratio;
    end
  else
    % �����������������Ϊ�Ӻ�Ѫ��
    Cap=C(i);
    JMat(FromNode+VesNum,FromNode+VesNum)=Cap*s;
    JMat(FromNode+VesNum,i)=1;
    JMat(FromNode+VesNum,inInd2)=-1;
  end
  
  % �������ڽڵ�
  outInd1=find(To(i)==From); % ��ѯ��������ö�Ѫ���յ��Ѫ��
  outInd2=find(To(i)==To);   % ��ѯ�����ö�Ѫ���յ��Ѫ��
  if length(outInd1)==2
    % ����ö�Ѫ����������Ѫ�ܣ�˵���ýڵ���һ���ֲ�ڵ�
    Cap=C(outInd1(1))+C(outInd1(2));
    JMat(ToNode+VesNum,ToNode+VesNum)=Cap*s;
    JMat(ToNode+VesNum,i)=-1;
    JMat(ToNode+VesNum,outInd1(1))=1;
    JMat(ToNode+VesNum,outInd1(2))=1;
  elseif length(outInd2)==2
    % �����2��Ѫ�ܵ�����յ㣬˵���ýڵ���һ����۽ڵ�
    outInd2(outInd2==i)=[];
    Cap=C(outInd1);
    JMat(ToNode+VesNum,ToNode+VesNum)=Cap*s;
    JMat(ToNode+VesNum,i)=-1;
    JMat(ToNode+VesNum,outInd2)=-1;
    JMat(ToNode+VesNum,outInd1)=1;
  elseif length(outInd1)==0
    % ������յ㲻���κ�Ѫ�ܵ���㣬˵���ýڵ���һ�������߽�
    % TODO(panqing): ֻ����������������õģ���������߽�Ӧͬ����߽�һ���Դ�
    if 1
      RL=(MeanP(i)-DeltaP(i)/2)*133/kg_ratio/(Flow(i)/60/1e12*m_ratio^2);
      switch DampFactor
        case DampVisc
          RL=RL*ViscRatio;
      end
      JMat(ToNode+VesNum,ToNode+VesNum)=1;
      JMat(ToNode+VesNum,outInd2)=-RL;
    else
      JMat(ToNode+VesNum,ToNode+VesNum)=1;
      JMat(ToNode+VesNum,outInd2)=-1;
    end
  else
    % �����������������Ϊ�Ӻ�Ѫ��
    Cap=C(outInd2);
    JMat(ToNode+VesNum,ToNode+VesNum)=Cap*s;
    JMat(ToNode+VesNum,i)=-1;
    JMat(ToNode+VesNum,outInd2)=1;
  end
end
% disp(cnt);
X=JMat\RHS;
