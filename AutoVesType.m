% 自动判断血管网络中血管类型 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1：动脉；2：毛细血管；3：静脉
function Type=AutoVesType(Boundary,From,To,Porder)
VelNum=length(From);   %血管总数
Type=zeros(VelNum,1);

%% 初步判断类型 %%%%
% 分叉为动脉
% 汇聚为静脉
% 其余为毛细血管
for i=1:length(Porder)
  j=Porder(i);
  ConvergeIndex=find(To==From(j));   %判断汇聚血管
  BifurIndex=find(From==To(j));   %判断分叉血管
  if length(BifurIndex)==2
    Type(j)=1;
  elseif (length(find(To==To(j)))==1 && length(BifurIndex)==1) %判断单连血管
    Type(j)=Type(BifurIndex);
  else
    Type(j)=2;
  end
  
  if length(ConvergeIndex)==2
    Type(j)=3;
  end
end

%% 边界血管修正 %%%%
for i=1:length(Boundary(:,1))
  bFromIndex=find(From==Boundary(i,1));
  if ~isempty(bFromIndex)   %动脉
    Type(bFromIndex)=1;
  end
end

%% 进一步修正 %%%%
for k=1:2  %循环次数，以减少计算顺序带来的差错
  for i=1:length(Porder)
    j=Porder(i);
    
    InIndex=find(To==From(j));
    OutIndex=find(From==To(j));
    InNum=length(InIndex);
    OutNum=length(OutIndex);
    if InNum==2   %2进
      if Type(InIndex(1))==3 || Type(InIndex(2))==3  %进血管中有静脉
        Type(j)=3;
      elseif Type(InIndex(1))==1 || Type(InIndex(2))==1  %进血管中有动脉
        Type(j)=2;
      else %进血管皆为毛细血管
        Type(j)=3;
      end
    elseif InNum==1 && OutNum==2 %1进2出
      if Type(OutIndex(1))==3 || Type(OutIndex(2))==3  %确保动静脉不直接相连
        Type(j)=2;
      else
        Type(j)=Type(InIndex);
      end
    elseif InNum==1 && OutNum==1 %1进1出
      if Type(InIndex)==3 %入血管为静脉
        Type(j)=3;
      else
        Type(j)=2;
      end
    elseif InNum==1 && OutNum==0 %1进0出
      if Type(InIndex)==1
        Type(j)=2;
      else
        Type(j)=3;
      end
    elseif OutNum==2 %0进2出
      if Type(OutIndex(1))==3 || Type(OutIndex(2))==3
        Type(j)=2;
      else
        Type(j)=1;
      end
    else %0进1出
      if Type(OutIndex(1))==3
        Type(j)=2;
      else
        Type(j)=1;
      end
    end
  end
end

end