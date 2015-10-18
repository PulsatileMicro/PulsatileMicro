%% 1D模型仿真结果分析
% Before running this script, what you have to do are:
% 1. Run main_1D.m (in the directory with main_1D.m file)
% 2. Enter the data directory (such as 546_qR_E1)
% 3. Run this script and click "add to path" but not "change directory"
clear;close all;clc;

%%%% 运行宏定义 %%%%
Macro   % 定义宏的文件

%%%% 是否检测脉搏波起点 %%%%
% 0:不检测 1:检测
PWV_flag=0;

%%%% 导入预处理的信息 %%%%
load ModelParam.mat;
load VesParam.mat;
load NetType.mat;
load Structure.mat;
load Order.mat;
load SegName.mat;
if NetTypeID<100
  load LinEquDatFile.mat;
  load FuncPara.mat;
end

global ADAP NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE N_PERIOD DT PERIOD ODEREL ODEABS IN_BOUND_TYPE RANDOM_PHASE LIN_PA VEL_PROFILE BG_LEN_RATIO BG_MASS_RATIO
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID SymNet_ID Junc_ID Single_ID Tree_ID
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP

PeriodNum=ModelParam(N_PERIOD);     % 仿真周期数
Period=ModelParam(PERIOD);          % 周期
dt=ModelParam(DT);                  % 步长
if ModelParam(ODESOLVER)==ONED_EXP
  dt=dt*10;                         % 对于ONED_EXP，步长为1e-4，但每隔10步输出一点
end
DispPeriodNum=2;                    % 最终用于分析，显示波形的周期数
NumHisPt=VesParam(20,1);            % 每段血管记录的点数
% NumHisPt=1;
VesNum=length(VesParam(1,:));       % 血管段数
colorOrder=64;                      % 颜色阶数
Len=VesParam(1,:)';                 % 血管长度
From=VesParam(10,:);                % 血管段起始点
To=VesParam(11,:);                  % 血管段终止点
len_ratio=ModelParam(BG_LEN_RATIO);
mass_ratio=ModelParam(BG_MASS_RATIO);
% t_plot: 用于显示的数据的下标
% all_plot: 所有数据的下标
t_plot=1/dt*Period*(PeriodNum-DispPeriodNum)+1:1/dt*Period*(PeriodNum);t_plot=round(t_plot);
all_plot=1/dt*Period*0+1:1/dt*Period*(PeriodNum);all_plot=round(all_plot);
% 读取his文件中的数据（均为SI单位制）
[MeanP1 MeanP2 MeanQ MeanU MeanA MeanVisc tAll PAll1 PAll2 UAll QAll AAll ViscAll]...
  =HemoAnal4Network(NumHisPt,NetTypeName,VesNum,t_plot,all_plot,ModelParam(ODESOLVER));
if ModelParam(ODESOLVER)==ONED_EXP || ModelParam(ODESOLVER)==ONED_IMP
  MeanP1=MeanP1/133;  % mmHg
  MeanP2=MeanP2/133;  % mmHg
  MeanQ=MeanQ*1e12*60;% nL/min
  MeanU=MeanU*1000;   % mm/s
  MeanA=MeanA*1e12;   % um^2
  PAll1=PAll1/133;    % mmHg
  PAll2=PAll2/133;    % mmHg
  QAll=QAll*1e12*60;  % nL/min
  UAll=UAll*1000;     % mm/s
  AAll=AAll*1e12;     % um^2
else
  QAll=QAll*60*1e3/len_ratio^3;
  MeanQ=MeanQ*60*1e3/len_ratio^3;
  PAll1=PAll1/133/mass_ratio;
  PAll2=PAll2/133/mass_ratio;
  MeanP1=MeanP1/133/mass_ratio;
  MeanP2=MeanP2/133/mass_ratio;
  for i=1:VesNum
    UAll(1,all_plot,i)=QAll(1,all_plot,i)./(0.25*pi.*Diam(i)^2)*1e6/len_ratio^3;
  end
end

%%%% 初始化各类分析参数变量 %%%%
PI_U=zeros(VesNum,NumHisPt);PI_P1=PI_U;PI_P2=PI_U;RPSI_U=PI_U;RPSI_P=PI_U;AMP_P=PI_U;AMP_U=PI_U;
OrgPTT_U=zeros(VesNum,1);
OrgPTT_P=zeros(VesNum,1);

% 每段血管的阻力和顺应性
Res=zeros(VesNum,1);
Cap=zeros(VesNum,1);

% 脉搏波起点
uInd=zeros(VesNum,2);

%% 计算 PI, Amplitude, PTT, RPSI %%%%
% PI & RPSI is computed based on the waveform in the middle point of each segment
% PTT is computed at the beginning of each segment
for j=1:VesNum
  % PI & Amplitude
  for i=1:NumHisPt
    % 所有参数均基于t_plot（稳定后）范围内的波形计算
    PI_U(j,i)=(max(UAll(i,t_plot,j))-min(UAll(i,t_plot,j)))./mean(UAll(i,t_plot,j));
    PI_P1(j,i)=(max(PAll1(i,t_plot,j))-min(PAll1(i,t_plot,j)))./mean(PAll1(i,t_plot,j));
    PI_P2(j,i)=(max(PAll2(i,t_plot,j))-min(PAll2(i,t_plot,j)))./mean(PAll2(i,t_plot,j));
    AMP_P(j,i)=max(PAll1(i,t_plot,j))-min(PAll1(i,t_plot,j));
    AMP_U(j,i)=max(UAll(i,t_plot,j))-min(UAll(i,t_plot,j));
  end

  if PWV_flag
    %   uInd(i,1)=Pulse_Start(PAll1(1,t_plot,i));
    %   uInd(i,2)=Pulse_Start(PAll1(3,t_plot,i));
    [signal,p_time,p_press,u_time,u_press]=findpu(PAll1(1,all_plot,i),1/dt,1);
    uInd(i,1)=u_time(end);
    [signal,p_time,p_press,u_time,u_press]=findpu(PAll1(3,all_plot,i),1/dt,1);
    uInd(i,2)=u_time(end);
  end
  % RPSI
  for i=1:NumHisPt
    % RPSI_U
    tmpU1=UAll(i,t_plot(1:20:end),j);
    tmpU=(tmpU1(1:end-2)+tmpU1(2:end-1)+tmpU1(3:end))/3;
    RPSI_U(j,i)=max(diff(tmpU)/20)/mean(tmpU)*1000;
%     tmpU=UAll(i,t_plot(1:end),j);
%     RPSI_U(j,i)=max(diff(tmpU)/20)*1000;
    % RPSI_P
    tmpP1=PAll1(i,t_plot(1:20:end),j);
    tmpP=(tmpP1(1:end-2)+tmpP1(2:end-1)+tmpP1(3:end))/3;
    RPSI_P(j,i)=max(diff(tmpP)/20)/mean(tmpP)*1000;
  end
end

% 生成变量的RGB分量，用于Graph软件绘制分布图
if ModelParam(ODESOLVER)==ONED_EXP || ModelParam(ODESOLVER)==ONED_IMP 
  HisInd=2;
else
  HisInd=1;
end
PI_U_RGB=ColorScal(PI_U(:,HisInd),colorOrder);
PI_P1_RGB=ColorScal(PI_P1(:,HisInd),colorOrder);
RPSI_U_RGB=ColorScal(RPSI_U(:,HisInd),colorOrder);
RPSI_P_RGB=ColorScal(RPSI_P(:,HisInd),colorOrder);
AMP_P_RGB=ColorScal(AMP_P(:,HisInd),colorOrder);
AMP_U_RGB=ColorScal(AMP_U(:,HisInd),colorOrder);
GenParaFile(PI_P1(:,HisInd),SegName,'PIP1');
GenParaFile(AMP_P(:,HisInd),SegName,'AMP_P');
GenParaFile(AMP_U(:,HisInd),SegName,'AMP_U');

%% 绘制变量的分布图
switch NetTypeID
  case {Net_546_ID,Net_546_Meas_ID}
    PlotTopol(PI_P1(:,HisInd),Net_546_ID,'PIP1');
    if ModelParam(ODESOLVER)==ONED_EXP || ModelParam(ODESOLVER)==ONED_IMP
      PlotTopol(PI_P2(:,HisInd),Net_546_ID,'PIP2');
    end
    PlotTopol(AMP_P(:,HisInd),Net_546_ID,'AMP_P');
    PlotTopol(AMP_U(:,HisInd),Net_546_ID,'AMP_U');
    InIndex=1;
    OutIndex=330;
    % 计算Typical pathway上的PWV, cm/s
    if PWV_flag
      PWV=zeros(4,1);
      PWV(1)=-sum(Len([1:6 52:53 75:76 82:84 88 91]))/(uInd(1,1)-uInd(91,2))/10;
      PWV(2)=-sum(Len([1:6]))/(uInd(1,1)-uInd(6,2))/10;
      PWV(3)=-sum(Len([52:53 75:76]))/(uInd(52,1)-uInd(76,2))/10;
      PWV(4)=-sum(Len([82:84 88 91]))/(uInd(82,1)-uInd(91,2))/10;
    end
  case Egg_818_ID
    
  case Net_389_ID
    PlotTopol(PI_P1(:,HisInd),Net_389_ID,'PIP');
    PlotTopol(PI_U(:,HisInd),Net_389_ID,'PIU');
    InIndex=100;
    OutIndex=87;
  case Net_913_ID
    PlotTopol(PI_P1(:,HisInd),Net_913_ID,'PIP');
    PlotTopol(PI_U(:,HisInd),Net_913_ID,'PIU');
    InIndex=376;
    OutIndex=395;
  case Egg_636_ID
    PlotTopol(PI_P1(:,HisInd),Egg_636_ID,'PIP');
    PlotTopol(PI_U(:,HisInd),Egg_636_ID,'PIU');
    InIndex=166;
    OutIndex=2;
  case SymNet_ID
    PIP_DampingRate=(PI_P1(1,1)-PI_P1(VesNum,3))/PI_P1(1,1);
    PIU_DampingRate=(PI_U(1,1)-PI_U(VesNum,3))/PI_U(1,1);
    figure;plot(PI_P1(OrderInd(:,1),2),'r.-');hold on;plot(PI_U(OrderInd(:,1),2),'o-');
%     figure;plot(RPSI_P(OrderInd(:,1),2),'r.-');hold on;plot(RPSI_U(OrderInd(:,1),2),'o-');
    figure;subplot(121);plot(PAll1(2,t_plot,1),'r-');subplot(122);plot(PAll1(2,t_plot,VesNum));
    PathPIP=PI_P1(OrderInd(:,1));PathPIU=PI_U(OrderInd(:,1));
  case Tree_ID
    if PWV_flag
      PWV=zeros(VesNum,1);
      Len=VesParam(1,:);
      for i=1:VesNum
        PWV(i)=(Len(i))/((uInd(i,2)-uInd(i,1))*dt); % cm/s
      end
    end
    DampingRate_U=PI_U(1,1)-PI_U(VesNum,3);
    DampingRate_P=(PI_P1(1,1)-PI_P1(VesNum,3))/PI_P1(1,1);
  case Junc_ID
    plot(PI_P1(:,2),'r.-');
    hold on;
    plot(PI_U(:,2),'o-');
end
PIP_DampingRate=(PI_P1(InIndex,HisInd)-PI_P1(OutIndex,HisInd))/PI_P1(InIndex,HisInd);
PIU_DampingRate=(PI_U(InIndex,HisInd)-PI_U(OutIndex,HisInd))/PI_U(InIndex,HisInd);
AMPP_DampingRate=(AMP_P(InIndex,HisInd)-AMP_P(OutIndex,HisInd))/AMP_P(InIndex,HisInd);
AMPU_DampingRate=(AMP_U(InIndex,HisInd)-AMP_U(OutIndex,HisInd))/AMP_U(InIndex,HisInd);

%% Reynolds Number & Womersley Number
ReyNum=1050*Diam*1e-6.*MeanU(:,HisInd)*1e-3./Visc(:,HisInd);
WomNum=Diam/2*1e-6.*sqrt(2*pi*75/60*1050./Visc(:,HisInd));

%% 计算每段血管的Order
if NetTypeID==Net_546_ID || NetTypeID==Net_546_Meas_ID
  [VesOrder OrderRange]=CalcVesOrder(From,To,830);
  % 绘制L/R Ratio关于Order的关系
  LR_Ratio=Len./Diam;
  for i=1:length(order_x)
    ind=(VesOrder==order_x(i));
    MeanLRRatio(i)=mean(LR_Ratio(ind));
  end
  for i=1:VesNum
    DelInd=(VesOrder==max(VesOrder));
    LR_Ratio(DelInd)=[];
    VesOrder(DelInd)=[];
  end
  %% Plot the Pressure,PI,PTT vs. Order
  % Init the plot vars
  MeanPIP_ord=[];MeanPIU_ord=[];MeanPTTP_ord=[];MeanPTTU_ord=[];MeanP_ord=[];
  StdPIP_ord=[];StdPIU_ord=[];StdPTTP_ord=[];StdPTTU_ord=[];StdP_ord=[];
  for i=1:length(order_x)
    ind=(AllOrder==order_x(i));
    MeanPIP_ord(i)=mean(PI_P1(ind));
    StdPIP_ord(i)=std(PI_P1(ind));
    MeanPIU_ord(i)=mean(PI_U(ind));
    StdPIU_ord(i)=std(PI_U(ind));
    MeanPTTP_ord(i)=mean(PTT_P(ind));
    StdPTTP_ord(i)=std(PTT_P(ind));
    MeanPTTU_ord(i)=mean(PTT_U(ind));
    StdPTTU_ord(i)=std(PTT_U(ind));
    MeanP_ord(i)=mean(MeanP(ind,2));
    StdP_ord(i)=std(MeanP(ind,2));
  end
  figure;hold on;whitebg('w');
  subplot(131)
  errorbar(order_x,MeanP_ord,StdP_ord,'k.-','LineWidth',2,'MarkerSize',24);
  xlabel('Generation number');
  ylabel('Pressure(mmHg)');
  subplot(132)
  errorbar(order_x,MeanPIP_ord,StdPIP_ord,'k.-','LineWidth',2,'MarkerSize',24);
  xlabel('Generation number');
  ylabel('PI_P1');
  legend('PI\_P');
  subplot(133);
  errorbar(order_x,MeanPTTP_ord,StdPTTP_ord,'k.-','LineWidth',2,'MarkerSize',24);
  xlabel('Generation number');
  ylabel('PTT_P(ms)');
  legend('PTT\_P');
end

if NetTypeID==Net_546_ID || NetTypeID==Net_546_Meas_ID
  %% 比较1D模型与稳态模型的仿真结果
  figure;
  subplot(221);hold on;
  plot(FuncPara(:,4),'r.-');
  plot(MeanP(:,2),'.-');
  ylabel('Pressure(mmHg)');
  legend('Pries data','1D data');
  subplot(222);hold on;
  plot(FuncPara(:,8),'r.-');
  plot(MeanQ(:,2),'.-');
  ylabel('Flow rate (nl/min)');
  legend('Pries data','1D data');
  subplot(223);hold on;
  plot(FuncPara(:,7),'r.-');
  % MeanU(MeanU(:,2)<0,2)=-MeanU(MeanU(:,2)<0,2);
  plot(MeanU(:,2),'.-');
  ylabel('Flow velocity (mm/s)');
  legend('Pries data','1D data');
  subplot(224);hold on;
  plot(FuncPara(:,3).*FuncPara(:,3)*pi/4,'r.-');
  plot(MeanA(:,2),'.-');
  ylabel('Area (um2)');
  legend('Pries data','1D data');
  
  ErrQ=sum(abs((MeanQ(:,2)-FuncPara(:,8))./FuncPara(:,8)))/VesNum;
  ErrP=sum(abs((MeanP(:,2)-FuncPara(:,4))./FuncPara(:,4)))/VesNum;
  ErrA=mean(abs(sqrt(4/pi*(MeanA(:,2))-FuncPara(:,3))./FuncPara(:,3)));
  ErrU=mean(abs((MeanU(:,2)-FuncPara(:,7))./FuncPara(:,7)));
end

%% 绘制每段血管的P，Q，U，A，Visc等波形
for j=1:VesNum
  h=figure;
  subplot(221);hold on;
  %     plot(tAll(1,t_plot,j), PAll1(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll1(2,t_plot,j), 'b', 'LineWidth', 2);
%   plot(tAll(2,round(OrgPTT_P(j)),j), PAll1(2,round(OrgPTT_P(j)),j), 'ro');
  %     plot(tAll(2,all_plot(t_plot(2):end),j), PpAll(t_plot(1):end,j), 'r', 'LineWidth', 2);
  %     plot(tAll(2,all_plot(t_plot(2):end),j), PnAll(t_plot(1):end,j), 'k', 'LineWidth', 2);
  %     plot(tAll(3,t_plot,j), PAll1(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Pressure(mmHg)');xlabel('t(s)');
  title(['Pressure Vessel ' num2str(SegName(j),'%04d')]);
  %     legend('StartPt', 'MidPt', 'EndPt');
  
  %     subplot(222);hold on;
  % %     plot(tAll(1,t_plot,j), QAll(1,t_plot,j), 'r', 'LineWidth', 2);
  %     plot(tAll(2,t_plot,j), QAll(2,t_plot,j), 'b', 'LineWidth', 2);
  % %     plot(tAll(3,t_plot,j), QAll(3,t_plot,j), 'k', 'LineWidth', 2);
  %     ylabel('Flow Rate(nl/min)');xlabel('t(s)');
  % %     legend('StartPt', 'MidPt', 'EndPt');
  %     title(['Flow Rate Vessel ' num2str(SegName(j),'%04d')]);
  
  subplot(222);hold on;
  %     plot(tAll(1,t_plot,j), UAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), UAll(2,t_plot,j), 'b', 'LineWidth', 2);
  %   plot(tAll(2,round(OrgPTT_U(j)),j), UAll(2,round(OrgPTT_U(j)),j), 'ro');
  %     plot(tAll(2,all_plot(t_plot(2):end),j), UpAll(t_plot(1):end,j), 'r', 'LineWidth', 2);
  %     plot(tAll(2,all_plot(t_plot(2):end),j), UnAll(t_plot(1):end,j), 'k', 'LineWidth', 2);
  %     plot(tAll(3,t_plot,j), UAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Velocity(mm/s)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Velocity Vessel ' num2str(SegName(j),'%04d')]);
  
  subplot(223);hold on;
  plot(tAll(2,t_plot,j), ViscAll(2,t_plot,j), 'b', 'LineWidth', 2);
  ylabel('Viscosity(mPas)');xlabel('t(s)');
  title(['Viscosity Vessel ' num2str(SegName(j),'%04d')]);
  
  subplot(224);hold on;
  %     plot(tAll(1,t_plot,j), PAll1(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll1(2,t_plot,j), 'b', 'LineWidth', 2);
  %     plot(tAll(3,t_plot,j), PAll1(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Area(um2)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Area Vessel ' num2str(SegName(j),'%04d')]);
  
%   saveas(h,['Vessel ' num2str(SegName(j),'%04d') '_' num2str(j,'%03d')]);
%   saveas(h,['Vessel ' num2str(SegName(j),'%04d') '_' num2str(j,'%03d')],'png');
  close(h);
end

%% Filter analysis
for i=1:VesNum
  for j=1:NumHisPt-1
    % 计算阻力
    % R=8uL/(pi*r^4)
    u=(ViscAll(j,t_plot,i)+ViscAll(j+1,t_plot,i))/2;
    L=Len(i)/1e6;
    r=(sqrt(AAll(j,t_plot,i)/pi)+sqrt(AAll(j+1,t_plot,i)/pi))/2/1e6;
    Res(i,:,j)=8*u.*L/pi./(r.^4);
    % 计算顺应性
    % C=AL/(rho*PWV^2)
    A=(AAll(j,t_plot,i)+AAll(j+1,t_plot,i))/2/1e12;
    PWV2=VesParam(4,i).*VesParam(3,i)./(2*1050*r);
    rh_ratio=r/VesParam(3,i);
    %     Cap(i,:,j)=A.*L./(1050*PWV2);
    % By Gross et al.
    Cap(i,:,j)=3*A.*(rh_ratio+1).^2*L/(VesParam(4,i)*(2*rh_ratio+1));
  end
end
f=1/(2*pi*Res.*Cap);
for i=1:VesNum
  for j=1:NumHisPt-1
    fmean(i,j)=mean(f(i,:,j));
    ResMean(i,j)=mean(Res(i,:,j));
    CapMean(i,j)=mean(Cap(i,:,j));
  end
end


