%% 1D模型仿真结果分析
% Before running this script, what you have to do are:
% 1. Run main_1D.m (in the directory with main_1D.m file)
% 2. Enter the data directory (such as 546_qR_E1)
% 3. Run this script and click "add to path" but not "change directory"
clear;close all;clc;

%%%% 运行宏定义 %%%%
Macro   % 定义宏的文件

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

PeriodNum=ModelParam(N_PERIOD);     % 仿真周期数
Period=ModelParam(PERIOD);          % 周期
dt=ModelParam(DT);                  % 步长
dt=dt*10;                           % 对于ONED_EXP，步长为1e-4，但每隔10步输出一点
DispPeriodNum=2;                    % 最终用于分析，显示波形的周期数
NumHisPt=VesParam(20,1);            % 每段血管记录的点数
VesNum=length(VesParam(1,:));       % 血管段数
colorOrder=64;                      % 颜色阶数
% t_plot: 用于显示的数据的下标
% all_plot: 所有数据的下标
t_plot=1/dt*Period*(PeriodNum-DispPeriodNum)+1:1/dt*Period*(PeriodNum);t_plot=round(t_plot);
all_plot=1/dt*Period*0+1:1/dt*Period*(PeriodNum);all_plot=round(all_plot);
% 读取his文件中的数据（均为SI单位制）
[MeanP MeanQ MeanU MeanA MeanVisc tAll PAll UAll QAll AAll ViscAll]...
  =HemoAnal4Network(NumHisPt,NetTypeName,VesNum,t_plot,all_plot);
MeanP=MeanP/133;    % mmHg
MeanQ=MeanQ*1e12*60;% nL/min
MeanU=MeanU*1000;   % mm/s
MeanA=MeanA*1e12;   % um^2
PAll=PAll/133;      % mmHg
QAll=QAll*1e12*60;  % nL/min
UAll=UAll*1000;     % mm/s
AAll=AAll*1e12;     % um^2

%%%% 初始化各类分析参数变量 %%%%
PI_U=zeros(VesNum,NumHisPt);PI_P=PI_U;RPSI_U=PI_U;RPSI_P=PI_U;Amp_P=PI_U;Amp_U=PI_U;
OrgPTT_U=zeros(VesNum,1);
OrgPTT_P=zeros(VesNum,1);

% 每段血管的阻力和顺应性
Res=zeros(VesNum,1);
Cap=zeros(VesNum,1);

%% %% 计算 PI, Amplitude, PTT, RPSI %%%%
% PI & RPSI is computed based on the waveform in the middle point of each segment
% PTT is computed at the beginning of each segment
for j=1:VesNum
  % PI & Amplitude
  for i=1:NumHisPt
    % 所有参数均基于t_plot（稳定后）范围内的波形计算
    PI_U(j,i)=(max(UAll(i,t_plot,j))-min(UAll(i,t_plot,j)))./mean(UAll(i,t_plot,j));
    PI_P(j,i)=(max(PAll(i,t_plot,j))-min(PAll(i,t_plot,j)))./mean(PAll(i,t_plot,j));
    Amp_P(j,i)=max(PAll(i,t_plot,j))-min(PAll(i,t_plot,j));
    Amp_U(j,i)=max(UAll(i,t_plot,j))-min(UAll(i,t_plot,j));
  end
  
  % PTT
  for i=1:NumHisPt
    % Foot-to-foot method
    if (mean(UAll(2,all_plot,j))>0 && j~=361)
%       [signal,p_time,p_press,u_time,u_press]=findpu(UAll(i,all_plot,j)+10,1000,1);
    else
%       [signal,p_time,p_press,u_time,u_press]=findpu(UAll(i,all_plot,j),1000,0);
    end
%     indu(i)=u_time(end);
%     [signal,p_time,p_press,u_time,u_press]=findpu(PAll(i,all_plot,j),1000,1);
%     [signal,p_time,p_press,u_time,u_press]=findpu(PAll(i,t_plot,j),1000,1);
%     close all;
%     indp(i)=u_time(end);
    %
    %     % Min value method
    % %     [maxValue indu(i)]=max(UAll(i,t_plot(1:round(end/Freq)),j));
    % %     [maxValue indp(i)]=max(PAll(i,t_plot(1:round(end/Freq)),j));
%     OrgPTT_P(j,i)=round(indp(i));
%     OrgPTT_U(j)=round(indu(i));
  end
  
  %   if indu(NumHisPt)<indu(1)
  %     indu(NumHisPt)=indu(NumHisPt)+round(10000/Freq);
  %   end
  %   if indp(NumHisPt)<indp(1)
  %     indp(NumHisPt)=indp(NumHisPt)+round(10000/Freq);
  %   end
  
  % RPSI
  for i=1:NumHisPt
    % RPSI_U
%     tmpU1=UAll(i,t_plot(1:20:end),j);
%     tmpU=(tmpU1(1:end-2)+tmpU1(2:end-1)+tmpU1(3:end))/3;
%     RPSI_U(j,i)=max(diff(tmpU)/20)/mean(tmpU)*1000;
    tmpU=UAll(i,t_plot(1:end),j);
    RPSI_U(j,i)=max(diff(tmpU)/20)*1000;
    % RPSI_P
    tmpP1=PAll(i,t_plot(1:20:end),j);
    tmpP=(tmpP1(1:end-2)+tmpP1(2:end-1)+tmpP1(3:end))/3;
    RPSI_P(j,i)=max(diff(tmpP)/20)/mean(tmpP)*1000;
  end
end

% 生成变量的RGB分量，用于绘制分布图
if NetTypeID==Net_546_ID
  PTT_P=OrgPTT_P-OrgPTT_P(1);
  PTT_U=OrgPTT_U-OrgPTT_U(1);
  tempInd=(PTT_P>=Period/dt);
  PTT_P(tempInd)=PTT_P(tempInd)-Period/dt;
  tempInd=(PTT_P<=-Period/dt*0.6);
  PTT_P(tempInd)=round(PTT_P(tempInd)+Period/dt);
end
% PTT_U_RGB=ColorScal(PTT_U,colorOrder);
% PTT_P_RGB=ColorScal(PTT_P(:,2),colorOrder);
PI_U_RGB=ColorScal(PI_U(:,2),colorOrder);
PI_P_RGB=ColorScal(PI_P(:,2),colorOrder);
RPSI_U_RGB=ColorScal(RPSI_U(:,2),colorOrder);
RPSI_P_RGB=ColorScal(RPSI_P(:,2),colorOrder);

switch NetTypeID
  case {Net_546_ID,Net_546_Meas_ID}
    %% 绘制变量的分布图
    % PlotMorph_546(PI_U);title('PI_U');
    PlotTopol(PI_P(:,2),1);
    % PlotMorph_546(PTT_U);title('PTT_U');
%     PlotMorph_546(PTT_P(:,2));title('PTT_P');
    % PlotMorph_546(RPSI_U);title('RPSI_U');
    % PlotMorph_546(RPSI_P);title('RPSI_P');
    PIP_DampingRate=(PI_P(1,1)-PI_P(330,2))/PI_P(1,1);
    PIU_DampingRate=(PI_U(1,1)-PI_U(330,2))/PI_U(1,1);
  case Egg_818_ID
    
  case Net_389_ID
    PlotTopol(PI_P(:,2),4);
    PlotTopol(PI_U(:,2),4);
    PIP_DampingRate=(PI_P(100,1)-PI_P(87,3))/PI_P(100,1);
    PIU_DampingRate=(PI_U(100,1)-PI_U(87,3))/PI_U(100,1);
  case Net_913_ID
    PIP_DampingRate=(PI_P(376,1)-PI_P(395,3))/PI_P(376,1);
    PIU_DampingRate=(PI_U(376,1)-PI_U(395,3))/PI_U(376,1);
    PlotTopol(PI_P(:,2),5);
  case Egg_636_ID
%     PlotTopol(MeanP,9);
    PlotTopol(PI_P(:,2),9);
    PIP_DampingRate=(PI_P(166,1)-PI_P(2,3))/PI_P(166,1);
    PIU_DampingRate=(PI_U(166,1)-PI_U(2,3))/PI_U(166,1);
  case SymNet_ID
    PIP_DampingRate=(PI_P(1,1)-PI_P(VesNum,3))/PI_P(1,1);
    PIU_DampingRate=(PI_U(1,1)-PI_U(VesNum,3))/PI_U(1,1);
    figure;plot(PI_P(OrderInd(:,1),2),'r.-');hold on;plot(PI_U(OrderInd(:,1),2),'o-');
%     figure;plot(RPSI_P(OrderInd(:,1),2),'r.-');hold on;plot(RPSI_U(OrderInd(:,1),2),'o-');
    figure;subplot(121);plot(PAll(2,t_plot,1),'r-');subplot(122);plot(PAll(2,t_plot,VesNum));
    PathPIP=PI_P(OrderInd(:,1));PathPIU=PI_U(OrderInd(:,1));
  case Tree_ID
    DampingRate_U=PI_U(1,1)-PI_U(VesNum,3);
    DampingRate_P=(PI_P(1,1)-PI_P(VesNum,3))/PI_P(1,1);
  case Junc_ID
    plot(PI_P(:,2),'r.-');
    hold on;
    plot(PI_U(:,2),'o-');
end

% 计算动脉树中的脉搏波传导速度
% 计算脉搏波起点
uInd=zeros(VesNum,2);
for i=1:VesNum
  uInd(i,1)=Pulse_Start(PAll(1,t_plot,i));
  uInd(i,2)=Pulse_Start(PAll(3,t_plot,i));
  %     if uInd(i,1)>800
  %       uInd(i,1)=uInd(i,1)-800;
  %     end
  %     if uInd(i,2)>800
  %       uInd(i,2)=uInd(i,2)-800;
  %     end
end
% 计算PWV
PWV=zeros(VesNum,1);
Len=VesParam(1,:);
for i=1:VesNum
  PWV(i)=(Len(i))/((uInd(i,2)-uInd(i,1))*dt); % cm/s
end

%% 计算546网络典型通路A的PWV cm/s
if NetTypeID==Net_546_ID || NetTypeID==Net_546_Meas_ID
  PWV(1)=-sum(Len([1:6 52:53 75:76 82:84 88 91]))/(uInd(1,1)-uInd(91,2))/10;
  PWV(2)=-sum(Len([1:6]))/(uInd(1,1)-uInd(6,2))/10;
  PWV(3)=-sum(Len([52:53 75:76]))/(uInd(52,1)-uInd(76,2))/10;
  PWV(4)=-sum(Len([82:84 88 91]))/(uInd(82,1)-uInd(91,2))/10;
end

%% 计算每段血管的Order
if NetTypeID==Net_546_ID
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
    MeanPIP_ord(i)=mean(PI_P(ind));
    StdPIP_ord(i)=std(PI_P(ind));
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
  ylabel('PI_P');
  legend('PI\_P');
  subplot(133);
  errorbar(order_x,MeanPTTP_ord,StdPTTP_ord,'k.-','LineWidth',2,'MarkerSize',24);
  xlabel('Generation number');
  ylabel('PTT_P(ms)');
  legend('PTT\_P');
end

if NetTypeID==Net_546_ID
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
  %     plot(tAll(1,t_plot,j), PAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll(2,t_plot,j), 'b', 'LineWidth', 2);
%   plot(tAll(2,round(OrgPTT_P(j)),j), PAll(2,round(OrgPTT_P(j)),j), 'ro');
  %     plot(tAll(2,all_plot(t_plot(2):end),j), PpAll(t_plot(1):end,j), 'r', 'LineWidth', 2);
  %     plot(tAll(2,all_plot(t_plot(2):end),j), PnAll(t_plot(1):end,j), 'k', 'LineWidth', 2);
  %     plot(tAll(3,t_plot,j), PAll(3,t_plot,j), 'k', 'LineWidth', 2);
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
  %     plot(tAll(1,t_plot,j), PAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll(2,t_plot,j), 'b', 'LineWidth', 2);
  %     plot(tAll(3,t_plot,j), PAll(3,t_plot,j), 'k', 'LineWidth', 2);
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

%% Reynolds Number & Womersley Number
% ReyNum=1050*Diam*1e-6.*MeanU(:,2)*1e-3./MeanVisc(:,2);
% WomNum=Diam/2*1e-6.*sqrt(2*pi*75/60*1050./MeanVisc(:,2));

% if NetTypeID==7
%   Result=[PI_P;PI_U;Amp_P;Amp_U];
% end