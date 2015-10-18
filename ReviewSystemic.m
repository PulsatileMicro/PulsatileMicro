%% 整体循环1D模型仿真结果分析
close all;clc;

PeriodNum=ModelParam(3);    % 仿真周期数
Period=ModelParam(2);
dt=ModelParam(1);
NumHisPt=3;
ODESolver=ModelParam(4);
DispPeriodNum=2;  % 最终用于分析，显示波形的周期数
% t_plot: 用于显示的数据的下标
% all_plot: 所有数据的下标
t_plot=1/dt*Period*(PeriodNum-DispPeriodNum)+1:1/dt*Period*(PeriodNum);t_plot=round(t_plot);
all_plot=1/dt*Period*0+1:1/dt*Period*(PeriodNum);all_plot=round(all_plot);
if ODESolver==ONED_EXP || ODESolver==ONED_IMP
  % 读取his文件中的数据（均为SI单位制）
  [MeanP MeanP2 MeanQ MeanU MeanA MeanVisc tAll PAll PAll2 UAll QAll AAll ViscAll]...
    =HemoAnal4Network(NumHisPt,NetTypeName,VesNum,t_plot,all_plot,ODESolver);
  MeanP=MeanP/133;    % mmHg
  MeanQ=MeanQ*1e6;    % mL/s
  MeanU=MeanU*1000;   % mm/s
  MeanA=MeanA*1e12;   % um^2
  PAll=PAll/133;
  QAll=QAll*1e6;
  UAll=UAll*1000;
  AAll=AAll*1e12;
else
  for j=1:VesNum
    fileName = [NetTypeName '_' int2str(j) '.his'];
    fid=fopen(fileName, 'r');
    C=textscan(fid, '%f %f %f %f', 'HeaderLines',2);
    Q=C{2};
    P1=C{3};  % start point
    P2=C{4};  % end point
    QAll(:,j)=Q;
    PAll1(:,j)=P1;
    PAll2(:,j)=P2;
    %     UAll(:,j)=QAll(:,j)./(0.25*pi*Diam(j).^2)*1e6/len_ratio^3;
  end
  PAll1=PAll1/133;
  PAll2=PAll2/133;
  QAll=QAll*1e6;
end

% 脉搏波起点分析
uInd=zeros(VesNum,1);
for i=1:VesNum
  if ODESolver==ONED_IMP || ODESolver==ONED_EXP
    uInd(i)=Pulse_Start(PAll(1,t_plot,i));
  else
    uInd(i)=Pulse_Start(PAll1(t_plot,i));
  end
end

% 计算cfPWV
% 从left common carotid(15)的base处，到right femoral artery(52)的起始位置
Dcf=sum(Len([14 18 27 28 35 37 39 41 43 50]));
cfPWV=Dcf./(uInd(52)-uInd(15))*1e3/1e2; % cm/s

% 绘制分析结果
figure;
if ODESolver==ONED_IMP || ODESolver==ONED_EXP
  plot(t_plot/1000,PAll(1,t_plot,52), 'LineWidth', 2);
  hold on;
  plot(t_plot/1000,PAll(1,t_plot,15), 'r', 'LineWidth', 2);
else
  plot(t_plot/1000,PAll1(t_plot,52), 'LineWidth', 2);
  hold on;
  plot(t_plot/1000,PAll1(t_plot,15), 'r', 'LineWidth', 2);
end
ylabel('Pressure(mmHg)');
xlabel('t(s)');
title('carotid vs. femoral pulse wave');
legend('femoral','carotid');
