% �Զ��ж�Ѫ��������Ѫ������ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1��������2��ëϸѪ�ܣ�3������
function Type=AutoVesType(Boundary,From,To,Porder)
VelNum=length(From);   %Ѫ������
Type=zeros(VelNum,1);

%% �����ж����� %%%%
% �ֲ�Ϊ����
% ���Ϊ����
% ����ΪëϸѪ��
for i=1:length(Porder)
  j=Porder(i);
  ConvergeIndex=find(To==From(j));   %�жϻ��Ѫ��
  BifurIndex=find(From==To(j));   %�жϷֲ�Ѫ��
  if length(BifurIndex)==2
    Type(j)=1;
  elseif (length(find(To==To(j)))==1 && length(BifurIndex)==1) %�жϵ���Ѫ��
    Type(j)=Type(BifurIndex);
  else
    Type(j)=2;
  end
  
  if length(ConvergeIndex)==2
    Type(j)=3;
  end
end

%% �߽�Ѫ������ %%%%
for i=1:length(Boundary(:,1))
  bFromIndex=find(From==Boundary(i,1));
  if ~isempty(bFromIndex)   %����
    Type(bFromIndex)=1;
  end
end

%% ��һ������ %%%%
for k=1:2  %ѭ���������Լ��ټ���˳������Ĳ��
  for i=1:length(Porder)
    j=Porder(i);
    
    InIndex=find(To==From(j));
    OutIndex=find(From==To(j));
    InNum=length(InIndex);
    OutNum=length(OutIndex);
    if InNum==2   %2��
      if Type(InIndex(1))==3 || Type(InIndex(2))==3  %��Ѫ�����о���
        Type(j)=3;
      elseif Type(InIndex(1))==1 || Type(InIndex(2))==1  %��Ѫ�����ж���
        Type(j)=2;
      else %��Ѫ�ܽ�ΪëϸѪ��
        Type(j)=3;
      end
    elseif InNum==1 && OutNum==2 %1��2��
      if Type(OutIndex(1))==3 || Type(OutIndex(2))==3  %ȷ����������ֱ������
        Type(j)=2;
      else
        Type(j)=Type(InIndex);
      end
    elseif InNum==1 && OutNum==1 %1��1��
      if Type(InIndex)==3 %��Ѫ��Ϊ����
        Type(j)=3;
      else
        Type(j)=2;
      end
    elseif InNum==1 && OutNum==0 %1��0��
      if Type(InIndex)==1
        Type(j)=2;
      else
        Type(j)=3;
      end
    elseif OutNum==2 %0��2��
      if Type(OutIndex(1))==3 || Type(OutIndex(2))==3
        Type(j)=2;
      else
        Type(j)=1;
      end
    else %0��1��
      if Type(OutIndex(1))==3
        Type(j)=2;
      else
        Type(j)=1;
      end
    end
  end
end

end