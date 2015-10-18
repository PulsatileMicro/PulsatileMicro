function [MeanP,Flow,Vel,DeltaP,Visc,Hd]=LinEqu(NetTypeID,DatMatrix,Boundary,DampPara)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit DampVisc

%% 1.数据解包
% 1.1 DatMatrix
% SegName=DatMatrix(:,1);
From=DatMatrix(:,2);
To=DatMatrix(:,3);
Len=DatMatrix(:,4);
Diam=DatMatrix(:,5);
% WallTh=DatMatrix(:,6);
% SegType=DatMatrix(:,7);
Visc=DatMatrix(:,8);
% E=DatMatrix(:,9);
VesNum=length(From);

% 1.2 Boundary
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% 1.3 DampPara
DampFactor=DampPara(1);
ViscRatio=DampPara(2);

%% 2.边界血管方向修正 %%%%
% 依据边界条件的符号修正
% 血流为正数->入边界->边界节点在From类中
% 血流为负数->出边界->边界节点在To类中
[From,To]=OriNodeModify(Boundary,From,To);

Qref=0.001; % TODO: Qref没用，这个程序需要重构
% Boundary(:,4)=0.45;
BHd=BoundaryInput(Boundary,From,To,Qref);

%% 3.单位调整
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s

%% 4. 迭代估算viscosity
% 最多Loop_Limit次循环
if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  Loop_Limit=1;
else
  Loop_Limit=50;
end
% 循环次数计数
Loop_Cnt=0;
% 初始Hd计算顺序. Porder, 正序. Norder, 逆序
Porder=1:VesNum;Norder=VesNum:1;
% 记录每次循环得到的Visc
DebugVisc=zeros(VesNum,Loop_Limit);
% 开始循环
while Loop_Cnt<Loop_Limit
  Loop_Cnt=Loop_Cnt+1;
  % 线性方程求解模块
  % 单位，输入：SI单位制，输出：P:mmHg, Q:nL/min
  [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam,Visc,Len,From,To);
  % 逆流处理模块
  [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju);
  % 调整Hd计算顺序
  [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,2,Eju);
  if Eju==1
    errorFlag=1;
    break;
  end
  % 计算Hd
  [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam,MeanFlowNew,0);
  
  % 自适应后的网络仿真不需要粘滞度调节
  if NetTypeID~=Net_546_ID && NetTypeID~=Net_122_ID
    %%粘滞度反馈
    umDiam=Diam.*1e6;   %um
    for i=1:VesNum
      Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
      switch DampFactor
        case DampVisc
          Visc(i)=Visc(i)*ViscRatio;
      end
    end
    % 循环进入第二次后(T>1)，开始校正Visc
    % 采用SOR方法校正，避免出现粘滞度振荡的情况
    if Loop_Cnt>1
      Visc=0.5*(Visc+DebugVisc(:,Loop_Cnt-1));
    end
  end
  
  % 记录每次迭代的仿真结果
  DebugVisc(:,Loop_Cnt)=Visc;
end

% 更新血管类型
if NetTypeID==Net_913_ID
  SegType=AutoVesType(Boundary,FromNew,ToNew,Porder);
elseif NetTypeID==Egg_818_ID
  ConstSegType=zeros(VesNum,1);
  % 部分血管段无法通过自动算法判断类型
  ConstSegType([216;544;547;559;561;47;57;70;80;143;681;683])=[3;ones(11,1)];
  SegType=autoVesTypeAdv(ConstSegType,FromNew,ToNew,Porder);
elseif NetTypeID==Egg_636_ID
  ConstSegType=zeros(VesNum,1);
  ConstSegType([167;168;143;144;560;561;562;212;580;581;583;410;615;449;514;452])=ones(16,1);
  ConstSegType([213;544;506;507;583;394;384;385;386;396])=2*ones(10,1);
  SegType=autoVesTypeAdv(ConstSegType,FromNew,ToNew,Porder);
  SegType(424)=3;
  SegType(428)=2;
end

%% 结果保存
Vel=MeanFlowNew/1e12/60./(pi*0.25*Diam.^2)*1e3;
Visc=Visc*1e3;
Flow=MeanFlowNew;
DeltaP=DeltaPNew;
% 保存steady state数据
save('LinEquDatFile.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
end
