%% 计算血管的Order
% Order如何定义
% Order只针对微动脉和毛细血管
function [VesOrder OrderRange]=CalcVesOrder(From,To,Inlet)
VesNum=length(From);
VesOrder=zeros(VesNum,1);
% 首先，找到主输入，将其Order定义为1
% VesOrder(From==Inlet)=1;
for i=1:VesNum
  tmpInd=i;
  order=1;
  while 1
    if From(tmpInd)==Inlet;
      VesOrder(i)=order;
      break;
    else
      inInd=find(From(tmpInd)==To);
      if length(inInd)==1
        tmpInd=inInd;
        order=order+1;
      elseif length(inInd)==2
        order=0;
        break;
      else
        order=0;
        break;
      end
    end
  end
end
VesOrder(VesOrder==0)=max(VesOrder)+1;
OrderRange=min(VesOrder):max(VesOrder);
end