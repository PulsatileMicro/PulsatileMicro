function [] = FREQ_ANAL_FUNC( input_args )
%%%% Ƶ������������� %%%%
% 0:���� 1:ѹ�� (TODO: ѹ�����벻��)
InputType=0;
% ������̬Ѫ������ѧ���� %%%%
% ���ж�̬����ǰ������Ҫ���ȷ�����̬Ѫ������ѧ����
% ����������Ϊ��̬ģ�����ñ߽���������Ҫ������߽�����
[SS_Press,SS_Flow,SS_Vel,SS_DeltaP,SS_Visc,SS_Hd]=...
  LinEqu(NetTypeID,DampFactor,VesNum,SegType,From,To,Diam,Len,Visc,Boundary,ViscRatio);
% ����Ѫ�������Ƶ������ %%%%
df=0.5;    % Ƶ�ʲ���
f=0:df:10;
[AllFlowX,AllPressX,NormAllFlowX,NormAllPressX,CutFreqFlow,CutFreqPress,AllDeltaPX,AllRatioPX]=...
  LinEquFreq(VesNum,From,To,Diam,Len,Visc,WallTh,E,Boundary,SS_Press,SS_Flow,SS_DeltaP,DampFactor,ViscRatio,InputType,f,df);
% ����ÿ��Ѫ�ܵ�Order
if NetTypeID==Net_546_ID || NetTypeID==Net_546_Meas_ID
  [AllOrder,order_x]=CalcVesOrder(From,To,830);
  MeanNormPressX=zeros(length(f),max(order_x));  % ��0Hz��10Hz��11��Ƶ��
  MeanNormFlowX=MeanNormPressX;StdNormPressX=MeanNormPressX;StdNormFlowX=MeanNormPressX;
  for j=1:length(f)
    for i=1:length(order_x)
      ind=find((AllOrder==order_x(i)));
      MeanNormPressX(j,i)=mean(NormAllPressX(j,ind));
      StdNormPressX(j,i)=std(NormAllPressX(j,ind));
      DelInd=find(NormAllFlowX(j,ind)>1 | NormAllFlowX(j,ind)<0);
      if ~isempty(DelInd)
        ind(DelInd)=[];
      end
      MeanNormFlowX(j,i)=mean(NormAllFlowX(j,ind));
      StdNormFlowX(j,i)=std(NormAllFlowX(j,ind));
    end
  end
end
MeanNormPressX=MeanNormPressX';
StdNormPressX=StdNormPressX';
MeanNormFlowX=MeanNormFlowX';
StdNormFlowX=StdNormFlowX';
end

