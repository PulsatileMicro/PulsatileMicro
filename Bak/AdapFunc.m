function y = AdapFunc(X,Porder,Norder,Qref,PO2ref,Tauref,Jo,Lref,t,Mo,...
  Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,WallTh,...
  From,To,BHd,BSO2in,BJm,BJc,DataArray)

% 边界数据
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

OrgDiam=Diam;
OrgVisc=Visc;

Loop2_Cnt=0;
kp=X(1);
km=X(2);
kc=X(3);
ks=X(4);
while Loop2_Cnt<Loop2_Adap_Num % 循环2，自适应因素计算
  Loop2_Cnt=Loop2_Cnt+1;
  Loop1_Cnt=0;
  Loop1_Visc_Num=50;
  
  % 数据记录矩阵
  DebugVisc=zeros(VesNum,Loop1_Visc_Num);
  
  while Loop1_Cnt<Loop1_Visc_Num   % 循环1，Visc反馈
    Loop1_Cnt=Loop1_Cnt+1;
    % 线性方程求解模块
    [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam/1e3,Visc,Len/1e3,From,To);
    % 逆流处理模块
    [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju);
    %   FromNew=From;ToNew=To;MeanFlowNew=MeanFlow;DeltaPNew=DeltaP;
    % 调整Hd计算顺序
    [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,2,Eju);
    if Eju==1
      errorFlag=1;
      break;
    end
    % 计算Hd
    [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam/1e3,MeanFlowNew,Eju);
    
    % 自适应后的网络仿真不需要粘滞度调节
    %       if NetTypeID~=Net_546_ID && NetTypeID~=Net_122_ID
    if 1
      %%粘滞度反馈
      umDiam=Diam.*1e3;   %um
      for i=1:VesNum
        Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
      end
      DebugVisc(:,Loop1_Cnt)=Visc;
      % 循环进入第二次后(Loop1_Cnt>1)，开始校正Visc
      % 采用SOR方法校正，避免出现粘滞度振荡的情况
      ViscModType=1;  %0：不修正；1：对Bifurcation修正；2：对All segments修正
      Alpha=0.5;
      DebugVisc=modifyViscosity(Porder,FromNew,ToNew,DebugVisc,Loop1_Cnt,ViscModType,Alpha,Eju);
      Visc=DebugVisc(:,Loop1_Cnt);
    end
    
    % 记录每次迭代的仿真结果
    %       DebugVisc(:,Loop1_Cnt)=Visc;
    %       DebugHd(:,Loop1_Cnt)=Hd;
    %       DebugP(:,Loop1_Cnt)=MeanP;
    %       DebugFlow(:,Loop1_Cnt)=MeanFlowNew;
    
    %% 循环结束模块 %%%%
    %MAE: max absolute error
    %LoopOutType:0-正常；1-收敛；2-达到步长未收敛；3-visc错误；4-线性方程错误
    [ViscMAE,LoopOutType,Loop1_Visc_Num]=ViscLoopTerminator(DebugVisc*1e3,Loop1_Visc_Num,Loop1_Cnt,Eju);
  end
  
  %PO2计算模块
  [PO2,SO2in,SO2mid,SO2out]=PO2Counter(Porder,BSO2in,FromNew,ToNew,Hd,Len,MeanFlowNew,Eju);
  %Sm计算模块
  [Sm,Jm]=SmCounter(Porder,BJm,FromNew,ToNew,Len,MeanFlowNew,PO2,Qref,PO2ref,t,Mo,Eju);
  %Sc计算模块
  [Sc,Jc]=ScCounter(Norder,BJc,FromNew,ToNew,Len,Sm,Jo,Lref,Eju);
  %Stau，Sp计算模块
  % TODO: 需要修改为Hypertension文献中的方程
  % [Stau,Sp,Tau]=SCounter(DeltaP,Diam,Len,MeanP,Tauref,Eju);
  %管径反馈模块
  % Stot=Stau+kp.*Sp+km.*(Sm+kc.*Sc)-ks;
  % Diam=Diam+Stot.*Diam.*t;
  
  % Hypertension文献方程
  Dm=Diam+WallTh/2;    %中壁管径(mm)
  Aw=WallTh.*2.*pi.*Dm; %中壁截面积(mm)
  w=WallTh*1e3;   %壁厚(um)
  tau=32*10^4/60.*Visc.*MeanFlowNew./(umDiam.^3);  %剪切力 dyn/cm2;
  %O=DeltaP*Diam/(2*WallTh)
  %以下计算公式中：DeltaP:mmHg; Diam:mm; WallTh:mm
  O=1333.2.*abs(DeltaP).*Diam./(2.*WallTh);    %周应力  dyn/cm2
  %参数 ref.Pries 2005
  ktd=1;
  kog=1;
  tauref=0.5598; %dyn/cm2
  Oref=32050;  %dyn/cm2
  wref=0.804;  %um
  kc=1.66;
  kmd=0.955;
  kmg=-0.374;
  ksd=3.077;
  ksg=0.0177;
  kwt=0.114;
  kwo=0.609;
  ee=0.1;
  %Rtd^2+Rod^2=1
  %Rtg^2+Rog^2=1
  Rtd=0.6;
  Rod=0.8;
  Rtg=0.6;
  Rog=0.8;
  Rw=0.2;
  Stm=ktd.*log10(tau/tauref+ee)./(1+kwt.*log10(w/wref+ee))+kmd.*(Sm+kc.*Sc)-ksd;
  Som=kog.*log10(O/Oref+ee)./(1+kwo.*log10(w/wref+ee))+kmg.*(Sm+kc.*Sc)-ksg;
  dDm=(Rtd.*Stm+Rod.*Som).*Dm.*t;
  dAw=Rw.*(Rtg.*Stm+Rog.*Som).*Aw.*t;
  newDm=Dm+dDm;     %mm
  newAw=Aw+dAw;     %mm2
  %迭代更新
  WallTh=newAw./(2.*pi.*newDm);   %mm
  Diam=newDm-WallTh/2;   %mm
  
  %% 主要参数记录
%   DebugStau(:,Loop2_Cnt)=Stau;
%   DebugSp(:,Loop2_Cnt)=Sp;
%   DebugSm(:,Loop2_Cnt)=Sm;
%   DebugSc(:,Loop2_Cnt)=Sc;
  DebugDiam(:,Loop2_Cnt)=Diam*1e3;
  DebugFlow(:,Loop2_Cnt)=MeanFlowNew;
  DebugP(:,Loop2_Cnt)=MeanP;
  DebugPO2(:,Loop2_Cnt)=PO2;
%   DebugTau(:,Loop2_Cnt)=Tau;
  DebugSO2in(:,Loop2_Cnt)=SO2in;
  DebugSO2out(:,Loop2_Cnt)=SO2out;
  DebugTHd(:,Loop2_Cnt)=Hd;
  
  %% 循环结束模块 %%%%
  %MAE: max absolute error
  %AdapLoopOutType:0-正常；1-收敛；2-达到步长未收敛；3-visc错误；4-线性方程错误
  %AccuracyType（精度控制）：0-最大绝对误差；1-最大相对误差；2-平均绝对误差
  AccuracyType=2;
  [DiamMAE,AdapLoopOutType,Loop1_Visc_Num]=AdapLoopTerminator(DebugDiam,Loop1_Visc_Num,Loop2_Cnt,Eju,AccuracyType);
  %异常错误输出参数修正
  if AdapLoopOutType>2
    if Loop2_Cnt==1
      Diam=OrgDiam;
    else  %错误跳出时，输出前一次循环数据，以便各分析参数的计算
      %为了计算Ev或Ed
      Diam=DebugDiam(:,Loop2_Cnt-1)./1e3;
      %为了计算加权的路径和能力耗散
      MeanFlow=abs(DebugFlow(:,Loop2_Cnt-1));
      % 为了计算氧需求
      SO2in=DebugSO2in(:,Loop2_Cnt-1);
      SO2out=DebugSO2out(:,Loop2_Cnt-1);
      Hd=DebugTHd(:,Loop2_Cnt-1);
      % 为了计算相对毛细血管压
      MeanP=DebugP(:,Loop2_Cnt-1);
    end
  end
  
  %% 输出显示 %%%%
  fprintf('Iterations:%3d  DiamMAE:%5.4e  viscIteration:%3d  ViscMAE:%5.4e.\n',Loop2_Cnt,DiamMAE,Loop1_Cnt,ViscMAE);
  if Loop1_Cnt==50
    break;
  end
end
%自适应仿真结果记录及分析%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 自适应灵敏度参数优化评价参数 %%%%
%Ed：管径变化差异；Ev：流速变化差异
umDiam=Diam*1e3;
%流速计算
Vel=4.*1e3.*MeanFlow./(60.*pi.*umDiam.*umDiam);   %血流速率mm/s
%Ref. Bettina 2009 (均方根误差）
Ev=sqrt(mean(4*(DataArray(:,8)-Vel).^2./(DataArray(:,8)+Vel).^2));

%对于非正常跳出的Ev惩罚
if AdapLoopOutType==2
  Ev=Ev*1.5;
else
  Ev=Ev*2;
end
y=Ev;
fprintf('Ev=%5.4e\n',y);
end


