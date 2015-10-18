% 红细胞比容计算 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BoundHd,FromNew,ToNew,Diam,MeanFlow,Eju)
%节点处红细胞总量守恒
%Hd边界初值设定
Hd=BoundHd;
VelNum=length(Porder);
FlowRatio=zeros(VelNum,1);
FQe=zeros(VelNum,1);

umDiam=Diam.*1e6;  %管径值单位转换，um
if Eju==0 %正常情况
  for i=1:length(Porder)
    j=Porder(i);
    ConvergeIndex=find(ToNew==FromNew(j));   %判断汇聚血管
    if length(ConvergeIndex)==2
      Hd(j)=(MeanFlow(ConvergeIndex(1))*Hd(ConvergeIndex(1))+MeanFlow(ConvergeIndex(2))*Hd(ConvergeIndex(2)))./(MeanFlow(j)); 
    end
    
    BifurIndex=find(FromNew==ToNew(j));   %判断分叉血管
    if length(BifurIndex)==2
      FlowRatio1=MeanFlow(BifurIndex(1))./(MeanFlow(j));
      FlowRatio2=MeanFlow(BifurIndex(2))./(MeanFlow(j));
      FQe1=PS_effect(umDiam(j),umDiam(BifurIndex(1)),umDiam(BifurIndex(2)),Hd(j),FlowRatio1);
      FQe2=PS_effect(umDiam(j),umDiam(BifurIndex(2)),umDiam(BifurIndex(1)),Hd(j),FlowRatio2);
      %特殊处理，由PhaseSeparation公式本身决定
      if ~(isreal(FQe1))==1 && ~(isreal(FQe2))==1
        if real(FQe1)>real(FQe2)
          FQe1=1;
          FQe2=0;
        else
          FQe1=0;
          FQe2=1;
        end
      end
      
      Hd(BifurIndex(1))=Hd(j)*FQe1/FlowRatio1;
      Hd(BifurIndex(2))=Hd(j)*FQe2/FlowRatio2;
%       Hd(BifurIndex(1))
%       Hd(BifurIndex(2))
      FlowRatio(BifurIndex(1))=FlowRatio1;
      FlowRatio(BifurIndex(2))=FlowRatio2;
      FQe(BifurIndex(1))=FQe1;
      FQe(BifurIndex(2))=FQe2; 
    elseif (length(find(ToNew==ToNew(j)))==1 && length(BifurIndex)==1) %判断单连血管   
      Hd(BifurIndex)=Hd(j);   
    else
      continue;
    end
  end
else
  Hd=zeros(VelNum,1);
end
end

