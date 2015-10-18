% 对流信号计算 上游->下游 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sm,Jmmid]=SmCounter(Porder,BJm,FromNeg,ToNeg,Len,MeanFlow,PO2,Qref,PO2ref,t,Mo,Eju)
%Jm=(Jm(in)+L*(1-PO2/PO2ref))*exp(-t/Mo);
%% 参数单位转换 %%%%
umLen=Len*1e3;  %单位转换um

%% 主要参数初值矩阵
VesNum=length(Porder);
Jm=1e-10*ones(VesNum,1);
Jmmid=zeros(VesNum,1);

%% 入边界各参数赋值 %%%%
Jmindex=find(BJm>0);  %光考BJm>0并不能判断所有边界，至少main feeding的Jm=0
Jm(1)=0;
for i=1:length(Jmindex)
  if PO2(Jmindex(i))<PO2ref
    Jm(Jmindex(i))=(BJm(Jmindex(i))+umLen(Jmindex(i))*(1-PO2(Jmindex(i))/PO2ref))*exp(-t/Mo);
    Jmmid(Jmindex(i))=BJm(Jmindex(i))+0.5*umLen(Jmindex(i))*(1-PO2(Jmindex(i))/PO2ref);
  else
    Jm(Jmindex(i))=BJm(Jmindex(i))*exp(-t/Mo);
    Jmmid(Jmindex(i))=BJm(Jmindex(i));
  end
end

if Eju==0  %正常情况
  for i=1:length(Porder)
    j=Porder(i);  %导入计算顺序矩阵
    ConvergeIndex=find(ToNeg==FromNeg(j));  %判断汇聚血管
    if length(ConvergeIndex)==2
      if PO2(j)<PO2ref
        Jm(j)=(Jm(ConvergeIndex(1))+Jm(ConvergeIndex(2))+umLen(j)*(1-PO2(j)/PO2ref))*exp(-t/Mo);
        Jmmid(j)=Jm(ConvergeIndex(1))+Jm(ConvergeIndex(2))+0.5*umLen(j)*(1-PO2(j)/PO2ref);
      else
        Jm(j)=(Jm(ConvergeIndex(1))+Jm(ConvergeIndex(2)))*exp(-t/Mo);
        Jmmid(j)=Jm(ConvergeIndex(1))+Jm(ConvergeIndex(2));
      end
    end
    BifurIndex=find(FromNeg==ToNeg(j));  %判断分叉血管
    if length(BifurIndex)==2
      if PO2(BifurIndex(1))>PO2ref&&PO2(BifurIndex(2))>PO2ref
        Jm(BifurIndex(1))=Jm(j)*MeanFlow(BifurIndex(1))/MeanFlow(j)*exp(-t/Mo);
        Jmmid(BifurIndex(1))=Jm(j)*MeanFlow(BifurIndex(1))/MeanFlow(j);
        Jm(BifurIndex(2))=Jm(j)*MeanFlow(BifurIndex(2))/MeanFlow(j)*exp(-t/Mo);
        Jmmid(BifurIndex(2))=Jm(j)*MeanFlow(BifurIndex(2))/MeanFlow(j);
      elseif PO2(BifurIndex(2))>PO2ref
        Jm(BifurIndex(1))=(Jm(j)*MeanFlow(BifurIndex(1))/MeanFlow(j)+umLen(BifurIndex(1))*(1-PO2(BifurIndex(1))/PO2ref))*exp(-t/Mo);
        Jmmid(BifurIndex(1))=Jm(j)*MeanFlow(BifurIndex(1))/MeanFlow(j)+0.5*umLen(BifurIndex(1))*(1-PO2(BifurIndex(1))/PO2ref);
        Jm(BifurIndex(2))=Jm(j)*MeanFlow(BifurIndex(2))/MeanFlow(j)*exp(-t/Mo);
        Jmmid(BifurIndex(2))=Jm(j)*MeanFlow(BifurIndex(2))/MeanFlow(j);
      elseif PO2(BifurIndex(1))>PO2ref
        Jm(BifurIndex(1))=Jm(j)*MeanFlow(BifurIndex(1))/MeanFlow(j)*exp(-t/Mo);
        Jmmid(BifurIndex(1))=Jm(j)*MeanFlow(BifurIndex(1))/MeanFlow(j);
        Jm(BifurIndex(2))=(Jm(j)*MeanFlow(BifurIndex(2))/MeanFlow(j)+umLen(BifurIndex(2))*(1-PO2(BifurIndex(2))/PO2ref))*exp(-t/Mo);
        Jmmid(BifurIndex(2))=Jm(j)*MeanFlow(BifurIndex(2))/MeanFlow(j)+0.5*umLen(BifurIndex(2))*(1-PO2(BifurIndex(2))/PO2ref);
      else
        Jm(BifurIndex(1))=(Jm(j)*MeanFlow(BifurIndex(1))/MeanFlow(j)+umLen(BifurIndex(1))*(1-PO2(BifurIndex(1))/PO2ref))*exp(-t/Mo);
        Jmmid(BifurIndex(1))=Jm(j)*MeanFlow(BifurIndex(1))/MeanFlow(j)+0.5*umLen(BifurIndex(1))*(1-PO2(BifurIndex(1))/PO2ref);
        Jm(BifurIndex(2))=(Jm(j)*MeanFlow(BifurIndex(2))/MeanFlow(j)+umLen(BifurIndex(2))*(1-PO2(BifurIndex(2))/PO2ref))*exp(-t/Mo);
        Jmmid(BifurIndex(2))=Jm(j)*MeanFlow(BifurIndex(2))/MeanFlow(j)+0.5*umLen(BifurIndex(2))*(1-PO2(BifurIndex(2))/PO2ref);
      end    
    elseif (length(find(ToNeg==ToNeg(j)))==1 && length(BifurIndex)==1)  %单连血管     
      if PO2(BifurIndex(1))<PO2ref
        Jm(BifurIndex(1))=(Jm(j)+umLen(BifurIndex(1))*(1-PO2(BifurIndex(1))/PO2ref))*exp(-t/Mo);
        Jmmid(BifurIndex(1))=Jm(j)+0.5*umLen(BifurIndex(1))*(1-PO2(BifurIndex(1))/PO2ref);
      else
        Jm(BifurIndex(1))=Jm(j)*exp(-t/Mo);
        Jmmid(BifurIndex(1))=Jm(j);
      end   
    else %分叉血管的子血管（已在分叉血管中被计算）
      continue;
    end
  end
else
  Jmmid=zeros(VesNum,1);
end
%计算对流信号
Sm=log10(1+Jmmid./(MeanFlow+Qref));
end
