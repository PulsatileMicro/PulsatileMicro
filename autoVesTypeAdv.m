% 自动判断血管网络中血管类型 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2013.11.11（改进）
%1：动脉；2：毛细血管；3：静脉
function y=autoVesTypeAdv(segType,From,To,Porder)
y=segType;
constTypeIndex=find(segType>0);  %人工确定的不变血管类型

%% 初步判断类型 %%%%
% 分叉为动脉
% 汇聚为静脉
% 其余为毛细血管
for i=1:length(Porder)
  j=Porder(i);
  if isempty(find(constTypeIndex==j,1))
    ConvergeIndex=find(To==From(j));   %判断汇聚血管
    BifurIndex=find(From==To(j));   %判断分叉血管
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

%% 进一步修正 %%%%
for k=1:2  %循环次数，以减少计算顺序带来的差错
  for i=1:length(Porder)
    j=Porder(i);
    if isempty(find(constTypeIndex==j,1))  %人工确定的血管不变
      InIndex=find(To==From(j));
      OutIndex=find(From==To(j));
      InNum=length(InIndex);
      OutNum=length(OutIndex);
      if InNum==2   %2进
        if y(InIndex(1))==3 || y(InIndex(2))==3  %进血管中有静脉
          y(j)=3;
        elseif y(InIndex(1))==1 || y(InIndex(2))==1  %进血管中有动脉
          if y(InIndex(1))==1 && y(InIndex(2))==1  %进血管都为动脉（特例）
            y(j)=1;
          else
            y(j)=2;
          end
        else %进血管皆为毛细血管
          y(j)=3;
        end
      elseif InNum==1 && OutNum==2 %1进2出
        if y(OutIndex(1))==3 || y(OutIndex(2))==3  %确保动静脉不直接相连
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
      elseif InNum==1 && OutNum==1 %1进1出
        if y(InIndex)==3 %入血管为静脉
          y(j)=3;
        else
          y(j)=2;
        end
      elseif InNum==1 && OutNum==0 %1进0出
        if y(InIndex)==1  %认为是动脉，毛细的情况需要人工判断
          y(j)=1;
        else
          y(j)=3;
        end
      elseif OutNum==2 %0进2出
        if y(OutIndex(1))==3 || y(OutIndex(2))==3
          if y(OutIndex(1))==3 && y(OutIndex(2))==3  %皆为静脉（特例）
            y(j)=3;
          else
            y(j)=3;
          end
        else
          y(j)=1;
        end
      else %0进1出
        if y(OutIndex(1))==3   %认为是静脉（毛细的情况需要人工判断）
          y(j)=3;
        else
          y(j)=1;
        end
      end
    end
  end
end

end