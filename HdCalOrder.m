% Ѫ���������˳��Ѱַ�����ں�����������Ĵ����źţ� %%%%%%%%%%%%%%%%%%%%%%%%%%
function [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNeg,ToNeg,P,N,DataType,Eju)
%������������ո�������˳�����Ѫ�ܣ�����һ��ѭ���м��������������������ֵ��
%�����д�����ѵļ���˳�򣬲�������������ʵ��Ѫ��˳��
%Porder������artery��vein�������Ѽ���˳��
%Norder������vein��artery�������Ѽ���˳��
VesNum=length(From);
Bond=Boundary(:,1);

% ����Porder
Pmarker=zeros(VesNum,1);
for i=1:length(Bond)
  FIndex= From==Bond(i);
  Pmarker(FIndex)=1;
end
% ��ֵ��Ϊ1��û��ֵ��Ϊ0.

if Eju==0 %�������
  count=0;
  Porder=zeros(VesNum,1);
  Pjudge=VesNum;
  while Pjudge>0
    Pjudge=Pjudge-1;
    for i=1:VesNum
      j=P(i);
      if j == 0
          continue;
      end
      ConvergeIndex=find(ToNeg==FromNeg(j));   %�жϻ��Ѫ��
      if length(ConvergeIndex)==2
        if Pmarker(ConvergeIndex(1))==0||Pmarker(ConvergeIndex(2))==0
          Pmarker(j)=0;
        else
          Pmarker(j)=1;
          %%Porder��¼�������Ѫ��˳��
          Index=find(Porder==ConvergeIndex(1), 1);
          if isempty(Index)
            count=count+1;
            Porder(count)=ConvergeIndex(1);
          end
          Index=find(Porder==ConvergeIndex(2), 1);
          if isempty(Index)
            count=count+1;
            Porder(count)=ConvergeIndex(2);
          end
          Index=find(Porder==j, 1);
          if isempty(Index)
            count=count+1;
            Porder(count)=j;
          end  
        end
      end
      
      BifurIndex=find(FromNeg==ToNeg(j));   %�жϷֲ�Ѫ��
      if length(BifurIndex)==2
        if Pmarker(j)==0
          Pmarker(j)=0;
        else
          Pmarker(BifurIndex(1))=1;
          Pmarker(BifurIndex(2))=1;
          %%Porder��¼�������Ѫ��˳��
          Index=find(Porder==j, 1);
          if isempty(Index)
            count=count+1;
            Porder(count)=j;
          end
          Index=find(Porder==BifurIndex(1), 1);
          if isempty(Index)
            count=count+1;
            Porder(count)=BifurIndex(1);
          end
          Index=find(Porder==BifurIndex(2), 1);
          if isempty(Index)
            count=count+1;
            Porder(count)=BifurIndex(2);
          end
        end
      elseif (length(find(ToNeg==ToNeg(j)))==1 && length(BifurIndex)==1) %�жϵ���Ѫ��
        if Pmarker(j)==0
          Pmarker(j)=0;
        else
          Pmarker(BifurIndex)=1;
          %%Porder��¼�������Ѫ��˳��
          Index=find(Porder==j, 1);
          if isempty(Index)
            count=count+1;
            Porder(count)=j;
          end
          Index=find(Porder==BifurIndex(1), 1);
          if isempty(Index)
            count=count+1;
            Porder(count)=BifurIndex(1);
          end
          
        end
      else
        continue;
      end
    end
    Index=find(Pmarker==0, 1);   %�ж�����
    if isempty(Index)==1
      Pjudge=0;
    else
      continue;
    end
  end
  
  if DataType<2  %��������Ӧ����
    %����Norder
    Nmarker=zeros(VesNum,1);
    for i=1:length(Bond)
      TIndex= To==Bond(i);
      Nmarker(TIndex)=1;
    end           %��ֵ��Ϊ1��û��ֵ��Ϊ0.
    
    From2=ToNeg;
    To2=FromNeg;
    count=0;
    Norder=zeros(VesNum,1);
    Njudge=VesNum;
    while Njudge>0
      Njudge=Njudge-1;
      for i=1:VesNum
        j=N(i);
        if j==0
            continue;
        end
        ConvergeIndex=find(To2==From2(j));   %�жϻ��Ѫ��
        if length(ConvergeIndex)==2         
          if Nmarker(ConvergeIndex(1))==0||Nmarker(ConvergeIndex(2))==0
            Nmarker(j)=0;
          else
            Nmarker(j)=1;
            %%Norder��¼�������Ѫ��˳��
            Index=find(Norder==ConvergeIndex(1), 1);
            if isempty(Index)
              count=count+1;
              Norder(count)=ConvergeIndex(1);
            end
            Index=find(Norder==ConvergeIndex(2), 1);
            if isempty(Index)
              count=count+1;
              Norder(count)=ConvergeIndex(2);
            end
            Index=find(Norder==j, 1);
            if isempty(Index)
              count=count+1;
              Norder(count)=j;
            end          
          end
        end      
        BifurIndex=find(From2==To2(j));   %�жϷֲ�Ѫ��
        if length(BifurIndex)==2     
          if Nmarker(j)==0
            Nmarker(j)=0;
          else
            Nmarker(BifurIndex(1))=1;
            Nmarker(BifurIndex(2))=1;
            %%Norder��¼�������Ѫ��˳��
            Index=find(Norder==j, 1);
            if isempty(Index)
              count=count+1;
              Norder(count)=j;
            end
            Index=find(Norder==BifurIndex(1), 1);
            if isempty(Index)
              count=count+1;
              Norder(count)=BifurIndex(1);
            end
            Index=find(Norder==BifurIndex(2), 1);
            if isempty(Index)
              count=count+1;
              Norder(count)=BifurIndex(2);
            end
          end       
        elseif (length(find(To2==To2(j)))==1 && length(BifurIndex)==1) %�жϵ���Ѫ��       
          if Nmarker(j)==0
            Nmarker(j)=0;
          else
            Nmarker(BifurIndex(1))=1;
            %%Norder��¼�������Ѫ��˳��
            Index=find(Norder==j, 1);
            if isempty(Index)
              count=count+1;
              Norder(count)=j;
            end
            Index=find(Norder==BifurIndex(1), 1);
            if isempty(Index)
              count=count+1;
              Norder(count)=BifurIndex(1);
            end          
          end
        else
          continue;
        end
      end    
      Index=find(Nmarker==0, 1);   %�ж�����
      if isempty(Index)==1
        Njudge=0;
      else
        continue;
      end
    end
  else
    Norder=N;  %ֻ��Ѫ������ѧ����
  end
else
  Porder=P;
  Norder=N;
end

if ~isempty(find(Porder==0, 1))||~isempty(find(Norder==0,1))
  Eju=1;
end

end