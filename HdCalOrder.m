% 血管网络计算顺序寻址（用于后续各沿网络的传输信号） %%%%%%%%%%%%%%%%%%%%%%%%%%
function [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNeg,ToNeg,P,N,DataType,Eju)
%输出参数（按照该数组编号顺序查找血管，可在一次循环中计算出整个网络的生理参数值）
%该序列代表最佳的计算顺序，并不代表网络真实的血流顺序
%Porder：正向，artery到vein方向的最佳计算顺序
%Norder：逆向，vein到artery方向的最佳计算顺序
VesNum=length(From);
Bond=Boundary(:,1);

% 正序Porder
Pmarker=zeros(VesNum,1);
for i=1:length(Bond)
  FIndex= From==Bond(i);
  Pmarker(FIndex)=1;
end
% 有值则为1；没有值则为0.

if Eju==0 %正常情况
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
      ConvergeIndex=find(ToNeg==FromNeg(j));   %判断汇聚血管
      if length(ConvergeIndex)==2
        if Pmarker(ConvergeIndex(1))==0||Pmarker(ConvergeIndex(2))==0
          Pmarker(j)=0;
        else
          Pmarker(j)=1;
          %%Porder记录被计算的血管顺序
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
      
      BifurIndex=find(FromNeg==ToNeg(j));   %判断分叉血管
      if length(BifurIndex)==2
        if Pmarker(j)==0
          Pmarker(j)=0;
        else
          Pmarker(BifurIndex(1))=1;
          Pmarker(BifurIndex(2))=1;
          %%Porder记录被计算的血管顺序
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
      elseif (length(find(ToNeg==ToNeg(j)))==1 && length(BifurIndex)==1) %判断单连血管
        if Pmarker(j)==0
          Pmarker(j)=0;
        else
          Pmarker(BifurIndex)=1;
          %%Porder记录被计算的血管顺序
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
    Index=find(Pmarker==0, 1);   %判断跳出
    if isempty(Index)==1
      Pjudge=0;
    else
      continue;
    end
  end
  
  if DataType<2  %含有自适应部分
    %逆序Norder
    Nmarker=zeros(VesNum,1);
    for i=1:length(Bond)
      TIndex= To==Bond(i);
      Nmarker(TIndex)=1;
    end           %有值则为1；没有值则为0.
    
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
        ConvergeIndex=find(To2==From2(j));   %判断汇聚血管
        if length(ConvergeIndex)==2         
          if Nmarker(ConvergeIndex(1))==0||Nmarker(ConvergeIndex(2))==0
            Nmarker(j)=0;
          else
            Nmarker(j)=1;
            %%Norder记录被计算的血管顺序
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
        BifurIndex=find(From2==To2(j));   %判断分叉血管
        if length(BifurIndex)==2     
          if Nmarker(j)==0
            Nmarker(j)=0;
          else
            Nmarker(BifurIndex(1))=1;
            Nmarker(BifurIndex(2))=1;
            %%Norder记录被计算的血管顺序
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
        elseif (length(find(To2==To2(j)))==1 && length(BifurIndex)==1) %判断单连血管       
          if Nmarker(j)==0
            Nmarker(j)=0;
          else
            Nmarker(BifurIndex(1))=1;
            %%Norder记录被计算的血管顺序
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
      Index=find(Nmarker==0, 1);   %判断跳出
      if isempty(Index)==1
        Njudge=0;
      else
        continue;
      end
    end
  else
    Norder=N;  %只有血流动力学部分
  end
else
  Porder=P;
  Norder=N;
end

if ~isempty(find(Porder==0, 1))||~isempty(find(Norder==0,1))
  Eju=1;
end

end