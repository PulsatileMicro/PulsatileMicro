function [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju)
%%输入参数
%From，To：血管节点关系
%MeanFlow：血流nl/min
%From1，To1为调整过的血管出入点矩阵，FlowMean保留原血流值，MeanFlow变为单位化后的血流的绝对值
%%输出参数
%From1，To1：逆流处调整后的血管节点关系
%FlowMean：血流nl/min
%MeanFlow：血流（绝对值）nl/min
FromNew=From;
ToNew=To;
if Eju==0
  FromNew=From;
  tempFrom=From;
  ToNew=To;
  tempTo=To;   %保持From，To的原始数据不变
  InvIndex=find(MeanFlow<0);  %判断血液逆流的血管编号
  for i=1:length(InvIndex)
    FromNew(InvIndex(i))=tempTo(InvIndex(i));
    ToNew(InvIndex(i))=tempFrom(InvIndex(i));
  end
else
  InvIndex=[];
end
MeanFlowNew=abs(MeanFlow);  %将流量取正，便于后续计算
DeltaPNew=abs(DeltaP);