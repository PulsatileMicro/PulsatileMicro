function [] = FREQ_ANAL_FUNC( input_args )
%%%% 频域分析输入类型 %%%%
% 0:流量 1:压力 (TODO: 压力输入不行)
InputType=0;
% 仿真稳态血流动力学特性 %%%%
% 进行动态仿真前，都需要首先仿真稳态血流动力学特性
% 仿真结果用于为动态模型设置边界条件，主要是输出边界条件
[SS_Press,SS_Flow,SS_Vel,SS_DeltaP,SS_Visc,SS_Hd]=...
  LinEqu(NetTypeID,DampFactor,VesNum,SegType,From,To,Diam,Len,Visc,Boundary,ViscRatio);
% 仿真血管网络的频响特性 %%%%
df=0.5;    % 频率步长
f=0:df:10;
[AllFlowX,AllPressX,NormAllFlowX,NormAllPressX,CutFreqFlow,CutFreqPress,AllDeltaPX,AllRatioPX]=...
  LinEquFreq(VesNum,From,To,Diam,Len,Visc,WallTh,E,Boundary,SS_Press,SS_Flow,SS_DeltaP,DampFactor,ViscRatio,InputType,f,df);
% 计算每段血管的Order
if NetTypeID==Net_546_ID || NetTypeID==Net_546_Meas_ID
  [AllOrder,order_x]=CalcVesOrder(From,To,830);
  MeanNormPressX=zeros(length(f),max(order_x));  % 从0Hz到10Hz，11个频率
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

