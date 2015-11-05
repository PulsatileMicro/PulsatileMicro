function [ObjFuncValue,AdapCoeff,OptDiam,OptWallTh,PO2,Sm,Sc] = AdapObjFunc(AdapCoeff,NetTypeID,AdapType,HdOrder,...
  AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray)
global WITH_WALL WITHOUT_WALL WITH_WALL_Cx WITHOUT_WALL_Cx WITH_WALL_Cx_PULSE

%% 1. 解包、预处理数据
% 1.1 解包血管网络数据
From=DatMatrix(:,2);
To=DatMatrix(:,3);
Len=DatMatrix(:,4);
Diam=DataArray(:,5);
WallTh=DatMatrix(:,6);
Visc=DatMatrix(:,8);
VesNum=length(From);
% NOTICE: 与其它子程序不同，此处单位调整成mm，参照王若帆的程序
Len=Len*1e-3;       % mm
Diam=Diam*1e-3;     % mm
Visc=Visc*1e-3;     % Pa.s
% WallTh=WallTh*1e-3; % mm

% 1.2 解包边界数据
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);
BHd=AdapBoundary(:,1);
BSO2in=AdapBoundary(:,2);
BJm=AdapBoundary(:,3);
BJc=AdapBoundary(:,4);

% 1.3 解包自适应仿真参数
dt=AdapPara(1);
M0=AdapPara(2);
Qref=AdapPara(3);
PO2ref=AdapPara(4);
J0=AdapPara(5);

% 1.4 Hd计算顺序
Porder=HdOrder(1,:);
Norder=HdOrder(2,:);

% 1.5 保存原始Diam和Visc
OrgDiam=Diam;

%% 2.设置自适应迭代参数
% 2.1 初始化迭代次数
Loop1_Visc_Limit=30;
Loop2_Adap_Limit=1000;

% 2.2 设置自适应参数
switch AdapType
  case WITH_WALL
    % 要优化的自适应系数
    kc=AdapCoeff(1);
    kmd=AdapCoeff(2);
    kmg=AdapCoeff(3);
    ksd=AdapCoeff(4);
    ksg=AdapCoeff(5);
    kwt=AdapCoeff(6);
    kwo=AdapCoeff(7);
    tauref=AdapCoeff(8);
    J0=AdapCoeff(9);
    Lref=AdapCoeff(10);
    Oref=AdapCoeff(11);
    wref=AdapCoeff(12);
  case WITHOUT_WALL
    % TODO
  case {WITH_WALL_Cx,WITH_WALL_Cx_PULSE}
    % 导入预存的自适应系数
    load AdapCoeff_Wall.mat;
    kc=AdapCoeff_Wall(1);
    kmd=AdapCoeff_Wall(2);
    kmg=AdapCoeff_Wall(3);
    ksd=AdapCoeff_Wall(4);
    ksg=AdapCoeff_Wall(5);
    kwt=AdapCoeff_Wall(6);
    kwo=AdapCoeff_Wall(7);
    tauref=AdapCoeff_Wall(8);
    J0=AdapCoeff_Wall(9);
    Lref=AdapCoeff_Wall(10);
    Oref=AdapCoeff_Wall(11);
    wref=AdapCoeff_Wall(12);
    % 要优化的Cx系数
    Cbasal=AdapCoeff(1);
    Cpeak=AdapCoeff(2);
    tau_peak=AdapCoeff(3);
    tau_width=AdapCoeff(4);
  case WITHOUT_WALL_Cx
    % TODO
end

%% 3. 开始自适应循环
% 初始化循环计数器
Loop2_Cnt=0;
while Loop2_Cnt<Loop2_Adap_Limit % 循环2，自适应循环迭代
  Loop2_Cnt=Loop2_Cnt+1;
  Loop1_Cnt=0;
  Loop1_Visc_Num=Loop1_Visc_Limit;
  % 数据记录矩阵
  DebugVisc=zeros(VesNum,Loop1_Visc_Num);
  FromNew=From;ToNew=To;
  while Loop1_Cnt<Loop1_Visc_Num   % 循环1，Visc迭代
    Loop1_Cnt=Loop1_Cnt+1;
    % 线性方程求解模块
    [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam/1e3,Visc,Len/1e3,FromNew,ToNew);
    % 逆流处理模块
    [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(FromNew,ToNew,MeanFlow,DeltaP,Eju);
    % 调整Hd计算顺序
    [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,1,Eju);
    if Eju==1
      errorFlag=1;
      break;
    end
    % 计算Hd
    [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam/1e3,MeanFlowNew,Eju);
    
    %%粘滞度反馈
    umDiam=Diam.*1e3;   % um
    for i=1:VesNum
      Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
    end
    DebugVisc(:,Loop1_Cnt)=Visc;
    % 循环进入第二次后(Loop1_Cnt>1)，开始校正Visc
    % 采用SOR方法校正，避免出现粘滞度振荡的情况
    ViscModType=1;  % 0：不修正；1：对Bifurcation修正；2：对All segments修正
    Alpha=0.5;
    DebugVisc=modifyViscosity(Porder,FromNew,ToNew,DebugVisc,Loop1_Cnt,ViscModType,Alpha,Eju);
    Visc=DebugVisc(:,Loop1_Cnt);
    
    %%%% 判断循环结束模块 %%%%
    % Loop1OutType:0-正常；1-收敛；2-达到步长未收敛；3-visc错误；4-线性方程错误
    [ViscMAE,Loop1OutType,Loop1_Visc_Num]=ViscLoopTerminator(DebugVisc*1e3,Loop1_Visc_Num,Loop1_Cnt,Eju);
  end
  
  % 血流动力学循环结束后，计算自适应因子
  switch AdapType
    case {WITH_WALL,WITH_WALL_Cx,WITH_WALL_Cx_PULSE}
      % Ref. Pries AR, Reglin B, and Secomb TW. Remodeling of blood vessels:
      % Responses of diameter and wall thickness to hemodynamic and metabolic stimuli.
      % Hypertension 46: 725-731, 2005.
      Dm=Diam+WallTh/2;       % 中壁管径(mm)
      Aw=WallTh.*2.*pi.*Dm;   % 中壁截面积(mm)
      w=WallTh*1e3;           % 壁厚(um)
      tau=32.*Visc.*MeanFlowNew/60*1e-12./(pi*Diam.^3)*1e9*10;  % 剪切力 dyne/cm2;
      %以下计算公式中：DeltaP:mmHg; Diam:mm; WallTh:mm
      O=133.*abs(MeanP).*Diam./(2.*WallTh)*10;                  % 周应力  dyne/cm2
      if AdapType==WITH_WALL_Cx || AdapType==WITH_WALL_Cx_PULSE
        Cex=Cx40exp(tau,Cbasal,Cpeak,tau_peak,tau_width);
        LrefVec=Lref*Cex.^0.5;
      else
        LrefVec=Lref.*ones(VesNum,1);
      end
      
      % 各部分自适应信号计算
      % PO2计算模块
      [PO2,SO2in,SO2mid,SO2out]=PO2Counter(Porder,BSO2in,FromNew,ToNew,Hd,Len,MeanFlowNew,Eju);
      % Sm计算模块
      [Sm,Jm]=SmCounter(Porder,BJm,FromNew,ToNew,Len,MeanFlowNew,PO2,Qref,PO2ref,dt,M0,Eju);
      % Sc计算模块
      [Sc,Jc]=ScCounter(Norder,BJc,FromNew,ToNew,Len,Sm,J0,LrefVec,Eju);
      
      ktd=1;
      kog=1;
      ee=1e-9;
      Rtd=0.8;Rod=0.05;Rtg=0.05;Rog=0.8;Rw=0.5;
      Stm=ktd.*log10(tau/tauref+ee)./(1+kwt.*log10(w/wref+ee))+kmd.*(Sm+kc.*Sc)-ksd;
      Som=kog.*log10(O/Oref+ee)./(1+kwo.*log10(w/wref+ee))+kmg.*(Sm+kc.*Sc)-ksg;
      dDm=(Rtd.*Stm+Rod.*Som).*Dm.*dt;
      dAw=Rw.*(Rtg.*Stm+Rog.*Som).*Aw.*dt;
      newDm=Dm+dDm;     % mm
      newAw=Aw+dAw;     % mm2
      % 迭代更新WallTh和Diam
      WallTh=newAw./(2.*pi.*newDm);   %mm
      Diam=newDm-WallTh/2;   %mm
      
      if AdapType==WITH_WALL_Cx_PULSE
        % 得到更新的WallTh和Diam后，生成动态模型仿真文件
        DatMatrix(:,2)=FromNew;
        DatMatrix(:,3)=ToNew;
        DatMatrix(:,5)=Diam*1e3;
        DatMatrix(:,6)=WallTh*1e3;
        DatMatrix(:,8)=Visc;
        SIM_PREP_FUNC(NetTypeID,Boundary,DatMatrix);
        cd('PulseAdapDIR');
        system('run.bat');
        % 运行完后进行分析
        cd('..');
      end
    case WITHOUT_WALL
      % TODO
    case WITHOUT_WALL_Cx
      % TODO
  end
  
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
  
  %% 自适应循环结束模块 %%%%
  switch AdapType
    case {WITH_WALL,WITH_WALL_Cx,WITH_WALL_Cx_PULSE}
      if mean(abs(dDm))<1e-5*dt && mean(abs(dAw))<1e-5*dt
        % 达到精度，跳出自适应循环
        AdapLoopOutType=1;
        break;
      elseif isnan(mean(abs(dDm))) || isnan(mean(abs(dAw)))
        AdapLoopOutType=4;
        break;
      elseif Loop1OutType==2
        break;
      end
      % 输出显示
      fprintf('AdapIt:=%d DiamMAE=%5.4e dDm:%5.4e dAw=%5.4e ViscIt:%3d.\n',Loop2_Cnt,1,mean(abs(dDm)),mean(abs(dAw)),Loop1_Cnt);
    case {WITHOUT_WALL,WITHOUT_WALL_Cx}
      % AdapLoopOutType:0-正常；1-收敛；2-达到步长未收敛；3-visc错误；4-线性方程错误
      % AccuracyType（精度控制）：0-最大管径绝对误差；1-最大管径相对误差；2-平均管径绝对误差；3-平均管径相对误差
      AccuracyType=3;
      [DiamMAE,AdapLoopOutType,Loop2_Adap_Num]=AdapLoopTerminator(DebugDiam,Loop2_Adap_Num,Loop2_Cnt,Eju,AccuracyType);
      % 异常错误输出参数修正
      if AdapLoopOutType>2
        if Loop2_Cnt==1
          Diam=OrgDiam;
        else  % 错误跳出时，输出前一次循环数据，以便各分析参数的计算
          % 为了计算Ev或Ed
          Diam=DebugDiam(:,Loop2_Cnt-1)./1e3;
        end
      end
      % 输出显示
      fprintf('AdapIt:=%d DiamMAE=%5.4e dDm:%5.4e dAw=%5.4e ViscIt:%3d.\n',Loop2_Cnt,1,mean(abs(dDm)),mean(abs(dAw)),Loop1_Cnt,ViscMAE);
  end
end

%% 4. 自适应灵敏度参数优化评价参数 %%%%
% Ed：管径变化差异；Ev：流速变化差异
umDiam=Diam*1e3;
% 流速计算
Vel=4.*1e3.*MeanFlowNew./(60.*pi.*umDiam.*umDiam);   %血流速率mm/s
if isempty(DataArray)
  Ed=sqrt(mean((OrgDiam-Diam).^2./(OrgDiam.^2)));
else
  Ev=sqrt(mean(4*(abs(DataArray(:,8))-Vel).^2./(abs(DataArray(:,8))+Vel).^2));
  Ed=sqrt(mean(4*(DataArray(:,5)-umDiam).^2./(DataArray(:,5).^2)));
end
% 对于非正常跳出的Ev惩罚
if Loop1OutType~=1 || Loop2_Cnt == 1000 || AdapLoopOutType~=1
  Ev=Ev*10;
  Ed=Ed*10;
end
% 记录优化目标函数值及相应的优化结果
ObjFuncValue=Ed;
OptDiam=Diam;
OptWallTh=WallTh;
fprintf('ObjFuncValue=%5.4e\n',ObjFuncValue);
