function [AdapPara,AllEv,AllWallTh,AllDiam,Visc,Sm,Sc]=AdapLinEqu(NetTypeID,DampFactor,VesNum,SegType,From,To,Diam,Len,Visc,WallTh,Boundary,AdapType,FuncPara,DataArray)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN E5 E02 Visc5 Visc02 DampE DampVisc

% 边界数据
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% 边界血管方向修正 %%%%
% 依据边界条件的符号修正
% 血流为正数->入边界->边界节点在From类中
% 血流为负数->出边界->边界节点在To类中
[From,To]=OriNodeModify(Boundary,From,To);

% 单位调整
% NOTICE: 与其它子程序不同，此处单位调整成mm，参照王若帆的程序
Len=Len*1e-3;   %mm
Diam=Diam*1e-3;   %mm
Visc=Visc*1e-3;   %Pa.s
WallTh=WallTh*1e-3;
% 保留原来的Diam，供计算Ed使用
OrgDiam=Diam;
OrgVisc=Visc;
OrgWallTh=WallTh;

% 自适应参数，初值
Qref=0.001;    %流量修正参数，对流模块
PO2ref=94.4;   %氧分压对照参数，对流模块
tauref=0.5598;  %剪切力修正参数，剪切力模块
Jo=6618;     %传导信号计算参数，传导模块
Lref=14292;    %血管长度衰减参数，传导模块
t=0.01;
Mo=1000;

if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  % if 0
  % 如果有自适应数据，那么只跑一次
  Loop1_ViscNum_Init=1;
  Loop2_Adap_Num=1;
  Loop3_PSO_Num=1;
  % 读取FuncPara里的数据作为边界
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo,FuncPara);
else
  % 最大循环次数 TODO(panqing):以收敛判断
  Loop1_ViscNum_Init=50;
  Loop2_Adap_Num_Init=2000;
  Loop3_PSO_Num=1;
  % 无FuncPara
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo);
end

% PSO参数
nPSO=1;
nParDim=7;
AdapPara=zeros(nParDim,nPSO,Loop3_PSO_Num);   % 所有粒子
v=zeros(nParDim,nPSO,Loop3_PSO_Num);          % 粒子速度
AllDiam=zeros(nPSO,VesNum);
AllWallTh=zeros(nPSO,VesNum);
InitAdapPara=[1.66,0.955,-0.374,3.077,0.0177,0.114,0.609];
for i=1:nPSO
%   AdapPara(:,i,1)=InitAdapPara+0.05*rand(1,nParDim)-0.05*rand(1,nParDim);
  AdapPara(:,i,1)=InitAdapPara;
end
c1=2; % 认知
c2=2; % 社会
r=1;  % 约束因子
ww=0.729;  % 惯性权重
AllEv=100*ones(nPSO,Loop3_PSO_Num);   % 所有Ev（目标函数），初值设为很大的100

% 记录每次循环得到的参数
DebugVisc=zeros(VesNum,Loop1_ViscNum_Init);
DebugHd=zeros(VesNum,Loop1_ViscNum_Init);
DebugP=zeros(VesNum,Loop1_ViscNum_Init);
DebugFlow=zeros(VesNum,Loop1_ViscNum_Init);

Loop3_Cnt=0;
% PSO迭代
while Loop3_Cnt<Loop3_PSO_Num
  Loop3_Cnt=Loop3_Cnt+1;
  Diam=OrgDiam;
  Visc=OrgVisc;
  WallTh=OrgWallTh;
  Loop2_Cnt=0;
  for j=1:nPSO
    kc=AdapPara(1,j,Loop3_Cnt);
    kmd=AdapPara(2,j,Loop3_Cnt);
    kmg=AdapPara(3,j,Loop3_Cnt);
    ksd=AdapPara(4,j,Loop3_Cnt);
    ksg=AdapPara(5,j,Loop3_Cnt);
    kwt=AdapPara(6,j,Loop3_Cnt);
    kwo=AdapPara(7,j,Loop3_Cnt);
      
    % 每次自适应迭代前，需将参数恢复为默认值
    Loop2_Cnt=0;
    Loop2_Adap_Num=Loop2_Adap_Num_Init;
    Diam=OrgDiam;
    Visc=OrgVisc;
    WallTh=OrgWallTh;
    while Loop2_Cnt<Loop2_Adap_Num % 循环2，自适应因素计算
      Loop2_Cnt=Loop2_Cnt+1;
      Loop1_Cnt=0;
      Loop1_Visc_Num=Loop1_ViscNum_Init;
      
      % 数据记录矩阵
      DebugConvergence=zeros(Loop1_Visc_Num-2,1);
      DebugVisc=zeros(VesNum,Loop1_Visc_Num);
      DebugHd=zeros(VesNum,Loop1_Visc_Num);  %记录循环2中的Hd变化
      DebugMeanFlow=zeros(VesNum,Loop1_Visc_Num);
      DebugDeltaP=zeros(VesNum,Loop1_Visc_Num);
      DebugMeanP=zeros(VesNum,Loop1_Visc_Num);
      DebugFQe=zeros(VesNum,Loop1_Visc_Num);
      
      LoopOutFlag=0;
      % 初始Hd计算顺序. Porder, 正序. Norder, 逆序
      Porder=1:VesNum;
      Norder=VesNum:-1:1;
%       Visc=OrgVisc;
      while Loop1_Cnt<Loop1_Visc_Num   % 循环1，Visc反馈
        Loop1_Cnt=Loop1_Cnt+1;
        % 线性方程求解模块
        [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam/1e3,Visc,Len/1e3,From,To);
        % 逆流处理模块
        [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju);
        % 调整Hd计算顺序
        [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,1,Eju);
        if Eju==1
          break;
        end
        % 计算Hd
        [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam/1e3,MeanFlowNew,Eju);
        
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
        
        %% 循环结束模块 %%%%
        %MAE: max absolute error
        %LoopOutType:0-正常；1-收敛；2-达到步长未收敛；3-visc错误；4-线性方程错误
        [ViscMAE,Loop1OutType,Loop1_Visc_Num]=ViscLoopTerminator(DebugVisc*1e3,Loop1_Visc_Num,Loop1_Cnt,Eju);
      end
      
      %PO2计算模块
      [PO2,SO2in,SO2mid,SO2out]=PO2Counter(Porder,BSO2in,FromNew,ToNew,Hd,Len,MeanFlowNew,Eju);
      %Sm计算模块
      [Sm,Jm]=SmCounter(Porder,BJm,FromNew,ToNew,Len,MeanFlowNew,PO2,Qref,PO2ref,t,Mo,Eju);
      %Sc计算模块
      [Sc,Jc]=ScCounter(Norder,BJc,FromNew,ToNew,Len,Sm,Jo,Lref,Eju);
      %Stau，Sp计算模块
      % TODO: 需要修改为Hypertension文献中的方程
      [Stau,Sp,Tau]=SCounter(DeltaP,Diam,Len,MeanP,tauref,Eju);
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
      Oref=32050;  %dyn/cm2
      wref=0.804;  %um
      ee=0.01;
      % 需满足：Rtd^2+Rod^2=1, Rtg^2+Rog^2=1
      Rtd=0.8;Rod=0;Rtg=0;Rog=0.6;Rw=0.2;
      Stm=ktd.*log10(tau/tauref+ee)./(1+kwt.*log10(w/wref+ee))+kmd.*(Sm+kc.*Sc)-ksd;
      Som=kog.*log10(O/Oref+ee)./(1+kwo.*log10(w/wref+ee))+kmg.*(Sm+kc.*Sc)-ksg;
      dDm=(Rtd.*Stm+Rod.*Som).*Dm.*t;
      dAw=Rw.*(Rtg.*Stm+Rog.*Som).*Aw.*t;
      newDm=Dm+dDm;     % mm
      newAw=Aw+dAw;     % mm2
      % 迭代更新WallTh和Diam
      WallTh=newAw./(2.*pi.*newDm);   %mm
      Diam=newDm-WallTh/2;   %mm
      
      %% 主要参数记录
      DebugStau(:,Loop2_Cnt)=Stau;
      DebugSp(:,Loop2_Cnt)=Sp;
      DebugSm(:,Loop2_Cnt)=Sm;
      DebugSc(:,Loop2_Cnt)=Sc;
      DebugDiam(:,Loop2_Cnt)=Diam*1e3;
      DebugFlow(:,Loop2_Cnt)=MeanFlowNew;
      DebugP(:,Loop2_Cnt)=MeanP;
      DebugPO2(:,Loop2_Cnt)=PO2;
      DebugTau(:,Loop2_Cnt)=Tau;
      DebugSO2in(:,Loop2_Cnt)=SO2in;
      DebugSO2out(:,Loop2_Cnt)=SO2out;
      DebugTHd(:,Loop2_Cnt)=Hd;
      
      %% 循环结束模块 %%%%
      %MAE: max absolute error
      %TLoopOutType:0-正常；1-收敛；2-达到步长未收敛；3-visc错误；4-线性方程错误
      %AccuracyType（精度控制）：0-最大绝对误差；1-最大相对误差；2-平均绝对误差；3-平均管径相对误差
      AccuracyType=3;
      [DiamMAE,TLoopOutType,Loop2_Adap_Num]=AdapLoopTerminator(DebugDiam,Loop2_Adap_Num,Loop2_Cnt,Eju,AccuracyType);
      %异常错误输出参数修正
      if TLoopOutType>2
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
      fprintf('PSOIter=%d, PSOParticle=%d, AdapIter=%d, DiamMAE:%5.4e,ViscIter=%d,ViscMAE:%5.4e.\n',Loop3_Cnt,j,Loop2_Cnt,DiamMAE,Loop1_Cnt,ViscMAE);
      if Loop1OutType~=1
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
    if isempty(DataArray)
      Ed=sqrt(mean((OrgDiam-Diam).^2./(OrgDiam.^2)));
    else
      Ev=sqrt(mean(4*(DataArray(:,8)-Vel).^2./(DataArray(:,8)+Vel).^2));
      Ed=sqrt(mean(4*(DataArray(:,5)-Diam*1e3).^2./(DataArray(:,5).^2)));
    end
    if LoopOutFlag==1
      %       Ev=Ev*10;
      Ed=Ed*10;
    end
    AllEv(j,Loop3_Cnt)=Ed;
    AllDiam(j,:)=Diam*1e3;
    AllWallTh(j,:)=WallTh*1e3;
    
    %     %% 自适应结果保存 %%%%
    %     AdaptationResult=[SegName DataArray(:,5) umDiam MeanP Tau FlowMean Hd PO2 Stot Stau Sp Sm Sc DeltaP];
    %     Parameters=[TLoopOutType AccuracyType DiamMAE Qref PO2ref Tauref Jo Lref kp km kc ks Ev Ed];
    %     AdaptData=[Parameters;AdaptationResult];
    %     save('AdaptData.mat','AdaptData');
    %     clear Vel umDiam;
  end
  % 计算适应度函数，更新下一次迭代的粒子
  % 1. 寻找历史搜索中的全局最优值
  [Value1,Ind1]=min(AllEv);
  [Value2,Ind2]=min(Value1);
  Pg=AdapPara(:,Ind1(Ind2),Ind2);
  % 2. 寻找历史搜索中，每个粒子的最优值
  for k=1:nPSO
    [Value,Ind]=min(AllEv(k,:));   % Ind为第k个粒子的历史最优值的序号
    Pi(:,k)=AdapPara(:,k,Ind);  % Pi(:,k)为第k个粒子的历史最优值
    v(:,k,Loop3_Cnt+1)=ww.*v(:,k,Loop3_Cnt)+c1*rand*(Pi(:,k)-AdapPara(:,k,Loop3_Cnt))+c2*rand*(Pg-AdapPara(:,k,Loop3_Cnt));
    AdapPara(:,k,Loop3_Cnt+1)=AdapPara(:,k,Loop3_Cnt)+r*v(:,k,Loop3_Cnt+1);
  end
end

