function [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju)
%%�������
%From��To��Ѫ�ܽڵ��ϵ
%MeanFlow��Ѫ��nl/min
%From1��To1Ϊ��������Ѫ�ܳ�������FlowMean����ԭѪ��ֵ��MeanFlow��Ϊ��λ�����Ѫ���ľ���ֵ
%%�������
%From1��To1���������������Ѫ�ܽڵ��ϵ
%FlowMean��Ѫ��nl/min
%MeanFlow��Ѫ��������ֵ��nl/min
FromNew=From;
ToNew=To;
if Eju==0
  FromNew=From;
  tempFrom=From;
  ToNew=To;
  tempTo=To;   %����From��To��ԭʼ���ݲ���
  InvIndex=find(MeanFlow<0);  %�ж�ѪҺ������Ѫ�ܱ��
  for i=1:length(InvIndex)
    FromNew(InvIndex(i))=tempTo(InvIndex(i));
    ToNew(InvIndex(i))=tempFrom(InvIndex(i));
  end
else
  InvIndex=[];
end
MeanFlowNew=abs(MeanFlow);  %������ȡ�������ں�������
DeltaPNew=abs(DeltaP);