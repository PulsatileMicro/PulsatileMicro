%% 0D模型仿真结果分析
% Before running this script, what you have to do are:
% 1. Enter the data directory (such as 546_qR_E1)
% 2. Run this script and click "add to path" but not "change directory"
clear;close all;clc;

%%%% 运行宏定义 %%%%
Macro   % 定义宏定义的文件

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
Period=ModelParam(PERIOD);
dt=ModelParam(DT);
len_ratio=ModelParam(BG_LEN_RATIO);
mass_ratio=ModelParam(BG_MASS_RATIO);
DispPeriodNum=2;                    % 最终用于分析，显示波形的周期数
NumHisPt=VesParam(20,1);
VesNum=length(VesParam(1,:));
colorOrder=64;                      % color map中颜色阶数
Len=VesParam(1,:)'/1e3/len_ratio;

% t_plot: 用于显示的数据的下标, all_plot: 所有数据的下标
t_plot=1/dt*Period*(PeriodNum-DispPeriodNum)+1:1/dt*Period*(PeriodNum);t_plot=round(t_plot);
all_plot=1/dt*Period*0+1:1/dt*Period*(PeriodNum);all_plot=round(all_plot);

% 读取his文件中的数据
QAll=zeros(length(all_plot),VesNum);UAll=QAll;PAll1=QAll;PAll2=QAll;
PI_P=zeros(VesNum,1);PI_U=PI_P;PIP1=PI_P;PIP2=PI_P;RPSI_P=PI_P;RPSI_U=PI_P;
MeanQ=PI_P;MeanU=PI_P;MeanPre=PI_P;P_Amp=PI_P;U_Amp=PI_P;MAP_P=PI_P;MAP_U=PI_P;
for j=1:VesNum
  if VesNum == 1
    fileName = [NetTypeName '.his'];
  else
    fileName = [NetTypeName '_' int2str(j) '.his'];
  end
  
  fid=fopen(fileName, 'r');
  C=textscan(fid, '%f %f %f %f', 'HeaderLines',2);
  Q=C{2};
  P1=C{3};  % start point
  P2=C{4};  % end point
  QAll(:,j)=Q;
  PAll1(:,j)=P1;
  PAll2(:,j)=P2;
  UAll(:,j)=QAll(:,j)./(0.25*pi*Diam(j).^2)*1e6/len_ratio^3;
  
  % PI
  PI_U(j)=(max(Q(t_plot))-min(Q(t_plot)))/mean(Q(t_plot));
  PIP1(j)=(max(P1(t_plot))-min(P1(t_plot)))/mean(P1(t_plot));
  PIP2(j)=(max(P2(t_plot))-min(P2(t_plot)))/mean(P2(t_plot));
  % 1/3Max+2/3Min
  MAP_P(j)=1/3*max(P1(t_plot))+2/3*min(P1(t_plot));
  MAP_U(j)=1/3*max(Q(t_plot))+2/3*min(Q(t_plot));
%   % RPSI
%   tmpU1=QAll(t_plot(1:20:end),j);
%   tmpU=(tmpU1(1:end-2)+tmpU1(2:end-1)+tmpU1(3:end))/3;
%   RPSI_U(j)=max(diff(tmpU)/20)/mean(tmpU)*1000;
%   tmpP1=PAll(t_plot(1:20:end),j);
%   tmpP=(tmpP1(1:end-2)+tmpP1(2:end-1)+tmpP1(3:end))/3;
%   RPSI_P(j)=max(diff(tmpP)/20)/mean(tmpP)*1000;
 
  % Mean parameters
  MeanQ(j)=mean(QAll(t_plot,j));
  MeanU(j)=mean(UAll(t_plot,j));
  if NetTypeID<100
    if abs(mean(P2(t_plot)))<0.1
      MeanPre(j)=(mean(P1(t_plot))+(MeanP(j)-DeltaP(j)/2)*133*mass_ratio)/2;
      PI_P(j)=PIP1(j);
    else
      MeanPre(j)=(mean(P1(t_plot))+mean(P2(t_plot)))/2;
      PI_P(j)=(PIP1(j)+PIP2(j))/2;
    end
  else
    MeanPre(j)=(mean(P1(t_plot))+mean(P2(t_plot)))/2;
    PI_P(j)=(PIP1(j)+PIP2(j))/2;
  end
  P_Amp(j)=PIP1(j)*MeanPre(j)/133; % mmHg
  U_Amp(j)=PI_U(j)*MeanU(j);       % mm/s
  
  % 绘制波形图
%   h=figure;
%   subplot(121);
%   plot(PAll1(t_plot,j));
%   subplot(122);
%   plot(QAll(t_plot,j));
%   saveas(h,['Vessel ' num2str(SegName(j),'%04d') '_' num2str(j,'%03d')],'png');
%   close(h);
  
  fclose(fid);
end
% MeanQ=MeanQ*1e6*60;
% MeanPre=MeanPre/1.33;
MeanQ=MeanQ*60*1e3/len_ratio^3;
MeanPre=MeanPre/133/mass_ratio;

% 计算动脉树中的脉搏波传导速度
% 计算脉搏波起点
uInd=zeros(VesNum,2);
for i=1:VesNum
%   uInd(i,1)=Pulse_Start(PAll1(t_plot,i));
%   uInd(i,2)=Pulse_Start(PAll2(t_plot,i));
  [signal,p_time,p_press,u_time,u_press]=findpu(PAll1(all_plot,i),1000,1);
  uInd(i,1)=u_time(end);
  [signal,p_time,p_press,u_time,u_press]=findpu(PAll2(all_plot,i),1000,1);
  uInd(i,2)=u_time(end);
end

switch NetTypeID
  case {Net_546_ID,Net_546_Meas_ID}
    % 衰减率
    PIP_DampingRate=(PIP1(1)-PIP1(330))/PIP1(1);
    PIU_DampingRate=(PI_U(1)-PI_U(330))/PI_U(1);
    % 衰减绝对值
    PIP_Damping=PIP1(1)-PIP1(330);
    PIU_Damping=PI_U(1)-PI_U(330);
    % 衰减率
    P_Amp_DampRate=(P_Amp(1)-P_Amp(330))/P_Amp(1);
    U_Amp_DampRate=(U_Amp(1)-U_Amp(330))/U_Amp(1);
    % 衰减量
    P_Amp_Damp=P_Amp(1)-P_Amp(330);
    U_Amp_Damp=U_Amp(1)-U_Amp(330);
    
    % Long Pathway
    ArtInd=[1 6 76 91];
    VenInd=[346 341 336 330];
    CapInd=219;
    x=[ArtInd CapInd VenInd];
    PathPI(:,1)=PI_P(x);
    PathPI(:,2)=PI_U(x);
    PathPI(:,3)=P_Amp(x);
    PathPI(:,4)=U_Amp(x);
    
    %% 计算546网络典型通路A的PWV cm/s
    if NetTypeID==Net_546_ID || NetTypeID==Net_546_Meas_ID
      PWV(1)=-sum(Len([1:6 52:53 75:76 82:84 88 91]))/(uInd(1,1)-uInd(91,2))/10*1e6;
      PWV(2)=-sum(Len([1:6]))/(uInd(1,1)-uInd(6,2))/10*1e6;
      PWV(3)=-sum(Len([52:53 75:76]))/(uInd(52,1)-uInd(76,2))/10*1e6;
      PWV(4)=-sum(Len([82:84 88 91]))/(uInd(82,1)-uInd(91,2))/10*1e6;
    end
    plot(t_plot,PAll1(t_plot,1));hold on;
    plot(uInd(1,1),PAll1(uInd(1,1),1),'ro');
    plot(t_plot,PAll2(t_plot,91),'--');hold on;
    plot(uInd(91,2),PAll2(uInd(91,2),91),'ko');
    
    figure;
    subplot(221);hold on;
%     plot(FuncPara(:,4),'r.-');
    plot(MeanP-min(MeanP),'r.-');
    plot(MeanPre-min(MeanPre),'.-');
    ylabel('Pressure(mmHg)');
    legend('Pries data','1D data');
    subplot(222);hold on;
%     plot(FuncPara(:,8),'r.-');
    plot(Flow,'r.-');
    plot(MeanQ,'.-');
    ylabel('Flow rate (nl/min)');
    legend('Pries data','1D data');
    subplot(223);hold on;
%     plot(FuncPara(:,7),'r.-');
    plot(Vel,'r.-');
    plot(MeanU,'.-');
    ylabel('Flow velocity (mm/s)');
    legend('Pries data','1D data');
  case Net_389_ID
    PIP_DampingRate=(PIP1(100)-PIP1(87))/PIP1(100);
    PIU_DampingRate=(PI_U(100)-PI_U(87))/PI_U(100);
    PIP_Damping=PIP1(100)-PIP1(87);
    PlotTopol(PI_P,NetTypeID);
  case Net_913_ID
    PIP_DampingRate=(PIP1(376)-PIP1(395))/PIP1(376);
    PIU_DampingRate=(PI_U(376)-PI_U(395))/PI_U(376);
    PIP_Damping=PIP1(376)-PIP1(395);
    PlotTopol(PI_P,NetTypeID);
  case Net_122_ID
    PIP_DampingRate=(PIP1(1)-PIP1(84))/PIP1(1);
    PIU_DampingRate=(PI_U(1)-PI_U(84))/PI_U(1);
    P_Amp_DampRate=(P_Amp(1)-P_Amp(84))/P_Amp(1);
    U_Amp_DampRate=(U_Amp(1)-U_Amp(84))/U_Amp(1);
    figure;
    subplot(221);hold on;
    plot(FuncPara(:,4)-min(FuncPara(:,4)),'r.-');
    plot(MeanPre-min(MeanPre),'.-');
    ylabel('Pressure(mmHg)');
    legend('Pries data','RC data');
    subplot(222);hold on;
    plot(FuncPara(:,8),'r.-');
    plot(MeanQ,'.-');
    ylabel('Flow rate (nl/min)');
    legend('Pries data','RC data');
    subplot(223);hold on;
    plot(FuncPara(:,7),'r.-');
    plot(MeanU,'.-');
    ylabel('Flow velocity (mm/s)');
    legend('Pries data','RC data');
  case Egg_CAM_ID
    PlotTopol(MeanP,Egg_CAM_ID);
    PlotTopol(PI_P,Egg_CAM_ID);
    PIP_DampingRate=(PIP1(974)-PIP1(6170))/PIP1(974);
    PIU_DampingRate=(PI_U(974)-PI_U(6170))/PI_U(974);
    PIP_DampingRate_Amp=(P_Amp(974)-P_Amp(6170))/P_Amp(974);
    PIU_DampingRate_Amp=(U_Amp(974)-U_Amp(6170))/U_Amp(974);
  case Sub_CAM_ID
    PlotTopol(PI_P,Sub_CAM_ID);
    PlotTopol(PI_U,Sub_CAM_ID);
    PIP_DampingRate=(PIP1(1)-PIP1(1241))/PIP1(1);
    PIU_DampingRate=(PI_U(1)-PI_U(1241))/PI_U(1);
  case Egg_818_ID
    PIP_DampingRate=(PIP1(761)-PIP1(357))/PIP1(761);
    PIU_DampingRate=(PI_U(761)-PI_U(357))/PI_U(761);
    PIP_DampingRate_Amp=(P_Amp(761)-P_Amp(357))/P_Amp(761);
    PIU_DampingRate_Amp=(U_Amp(761)-U_Amp(357))/U_Amp(761);
    PlotTopol(PI_P,Egg_818_ID);
    PlotTopol(PI_U,Egg_818_ID);
  case SymNet_ID
    SymInd=zeros(order,1);
    for i=1:order*2-1
      if i<=order
        SymInd(i)=2^(i-1);
      else
        SymInd(i)=SymInd(i-1)+SymInd(2*order-i+1);
      end
    end
    figure;plot(PIP1(SymInd),'r.-');hold on;plot(PI_U(SymInd),'o-');ylim([0 1.2]);legend('PI_P','PI_U');
%     figure;plot(P_Amp(SymInd),'r.-');hold on;plot(U_Amp(SymInd),'o-');
%     figure;plot(MeanPre(SymInd),'r.-');hold on;plot(MeanQ(SymInd),'o-');
    figure;subplot(121);plot(PAll1(t_plot,SymInd(1)));subplot(122);plot(PAll1(t_plot,SymInd(end)));
%     PI_All=[PIP1(SymInd(1)),PIP1(SymInd(end));PI_U(SymInd(1)),PI_U(SymInd(end))]';
    PI_All=[PIP1(SymInd),PI_U(SymInd),P_Amp(SymInd),U_Amp(SymInd),...
      P_Amp(SymInd)./P_Amp(1),U_Amp(SymInd)./U_Amp(1),MeanPre(SymInd)./MeanPre(1),MeanQ(SymInd)./MeanQ(1)];
    PIP_DampingRate=(PIP1(1)-PIP1(SymInd(end)))/PIP1(1);
    PIU_DampingRate=(PI_U(1)-PI_U(SymInd(end)))/PI_U(1);
    PIP_Damping=PIP1(1)-PIP1(SymInd(end));
  case Tree_ID
    SymInd=zeros(order,1);
    for i=1:order
      SymInd(i)=2^(i-1);
    end
    figure;plot(PIP1(SymInd),'r.-');hold on;plot(PI_U(SymInd),'o-');ylim([0 1.2]);
    figure;subplot(121);plot(PAll1(t_plot,SymInd(1)));subplot(122);plot(PAll1(t_plot,SymInd(end)));
    PI_All=[PIP1(SymInd),PI_U(SymInd),P_Amp(SymInd),U_Amp(SymInd)];
    PIP_DampingRate=(PIP1(1)-PIP1(SymInd(end)))/PIP1(1);
    PIU_DampingRate=(PI_U(1)-PI_U(SymInd(end)))/PI_U(1);
  case Junc_ID
    figure;plot(PIP1,'r.-');hold on;plot(PI_U,'o-');
%     ylim([0 1.2]);
    PIP_DampingRate=(PIP1(1)-PIP1(VesNum))/PIP1(1);
    PIU_DampingRate=(PI_U(1)-PI_U(VesNum))/PI_U(1);
  case Single_ID
end

%% 结果显示
fprintf(1,'仿真参数设置\n');
fprintf(1,'PIP Damping Rate: %.2f%%\n',PIP_DampingRate*100);
fprintf(1,'PIU Damping Rate: %.2f%%\n',PIU_DampingRate*100);
fprintf(1,'PIP Damping: %.2f\n',PIP_Damping);
fprintf(1,'PIU Damping: %.2f\n',PIU_Damping);
fprintf(1,'PAmp Damping Rate: %.2f%%\n',P_Amp_DampRate*100);
fprintf(1,'UAmp Damping Rate: %.2f%%\n',U_Amp_DampRate*100);
