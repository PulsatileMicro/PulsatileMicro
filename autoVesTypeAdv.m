% �Զ��ж�Ѫ��������Ѫ������ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2013.11.11���Ľ���
%1��������2��ëϸѪ�ܣ�3������
function y=autoVesTypeAdv(segType,From,To,Porder)
y=segType;
constTypeIndex=find(segType>0);  %�˹�ȷ���Ĳ���Ѫ������

%% �����ж����� %%%%
% �ֲ�Ϊ����
% ���Ϊ����
% ����ΪëϸѪ��
for i=1:length(Porder)
  j=Porder(i);
  if isempty(find(constTypeIndex==j,1))
    ConvergeIndex=find(To==From(j));   %�жϻ��Ѫ��
    BifurIndex=find(From==To(j));   %�жϷֲ�Ѫ��
    if length(BifurIndex)==2
      y(j)=1;
    else
      y(j)=2;
    end
    
    if length(ConvergeIndex)==2
      y(j)=3;
    end
  end
end

%% ��һ������ %%%%
for k=1:2  %ѭ���������Լ��ټ���˳������Ĳ��
  for i=1:length(Porder)
    j=Porder(i);
    if isempty(find(constTypeIndex==j,1))  %�˹�ȷ����Ѫ�ܲ���
      InIndex=find(To==From(j));
      OutIndex=find(From==To(j));
      InNum=length(InIndex);
      OutNum=length(OutIndex);
      if InNum==2   %2��
        if y(InIndex(1))==3 || y(InIndex(2))==3  %��Ѫ�����о���
          y(j)=3;
        elseif y(InIndex(1))==1 || y(InIndex(2))==1  %��Ѫ�����ж���
          if y(InIndex(1))==1 && y(InIndex(2))==1  %��Ѫ�ܶ�Ϊ������������
            y(j)=1;
          else
            y(j)=2;
          end
        else %��Ѫ�ܽ�ΪëϸѪ��
          y(j)=3;
        end
      elseif InNum==1 && OutNum==2 %1��2��
        if y(OutIndex(1))==3 || y(OutIndex(2))==3  %ȷ����������ֱ������
          if y(InIndex)==1
            y(j)=2;
          else
            y(j)=3;
          end
        else
          y(j)=y(InIndex);
          y(OutIndex(1))=y(InIndex);
          y(OutIndex(2))=y(InIndex);
        end
      elseif InNum==1 && OutNum==1 %1��1��
        if y(InIndex)==3 %��Ѫ��Ϊ����
          y(j)=3;
        else
          y(j)=2;
        end
      elseif InNum==1 && OutNum==0 %1��0��
        if y(InIndex)==1  %��Ϊ�Ƕ�����ëϸ�������Ҫ�˹��ж�
          y(j)=1;
        else
          y(j)=3;
        end
      elseif OutNum==2 %0��2��
        if y(OutIndex(1))==3 || y(OutIndex(2))==3
          if y(OutIndex(1))==3 && y(OutIndex(2))==3  %��Ϊ������������
            y(j)=3;
          else
            y(j)=3;
          end
        else
          y(j)=1;
        end
      else %0��1��
        if y(OutIndex(1))==3   %��Ϊ�Ǿ�����ëϸ�������Ҫ�˹��жϣ�
          y(j)=3;
        else
          y(j)=1;
        end
      end
    end
  end
end

end