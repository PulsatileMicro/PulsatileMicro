%% 单根血管结果分析
close all;clc;
% 从his文件中读取血流动力学仿真结果
PeriodNum=6;    % 仿真周期数
t_plot=1+800*(PeriodNum-1):800*PeriodNum;
all_plot=1:800*PeriodNum;
NumHisPt=VesParam(21,1);
[MeanP MeanQ MeanU MeanA MeanVisc tAll PAll UAll QAll AAll ViscAll]=HemoAnal4Network(NumHisPt,VesType,VesNum,t_plot,all_plot);
MeanP=MeanP/133;MeanQ=MeanQ*1e12*60;MeanU=MeanU*1000;MeanA=MeanA*1e12;PAll=PAll/133;QAll=QAll*1e12*60;UAll=UAll*1000;AAll=AAll*1e12;

% 绘制每段血管的P, Q, U, A曲线
for j=1:VesNum
  h=figure;
  subplot(221);hold on;
%   plot(tAll(1,t_plot,j), PAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll(2,t_plot,j), 'b', 'LineWidth', 2);
%   plot(tAll(3,t_plot,j), PAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Pressure(mmHg)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Pressure Vessel ' num2str(SegName(j),'%04d')]);
  %     ylim([0 3])
  
  subplot(222);hold on;
  plot(tAll(1,t_plot,j), QAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), QAll(2,t_plot,j), 'b', 'LineWidth', 2);
  plot(tAll(3,t_plot,j), QAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Flow Rate(nl/min)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Flow Rate Vessel ' num2str(SegName(j),'%04d')]);
  
  subplot(223);hold on;
  plot(tAll(1,t_plot,j), UAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), UAll(2,t_plot,j), 'b', 'LineWidth', 2);
  plot(tAll(3,t_plot,j), UAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Velocity(mm/s)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Velocity Vessel ' num2str(SegName(j),'%04d')]);
  %     ylim([0 10])
  
  subplot(224);hold on;
  plot(tAll(1,t_plot,j), PAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll(2,t_plot,j), 'b', 'LineWidth', 2);
  plot(tAll(3,t_plot,j), PAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Area(um2)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Area Vessel ' num2str(SegName(j),'%04d')]);
  
%   print(h,'-dpng',['Vessel ' num2str(SegName(j),'%04d')]);
%   close(h);
end

%% 计算PI和RPSI
PI_U=zeros(NumHisPt,1);PI_P=zeros(NumHisPt,1);PP_P=zeros(NumHisPt,1);PP_U=zeros(NumHisPt,1);RPSI_U=zeros(NumHisPt,1);RPSI_P=zeros(NumHisPt,1);
for i=1:NumHisPt
  PI_U(i)=(max(UAll(i,t_plot))-min(UAll(i,t_plot)))./mean(UAll(i,t_plot));
  PI_P(i)=(max(PAll(i,t_plot))-min(PAll(i,t_plot)))./mean(PAll(i,t_plot));
  PP_P(i)=max(UAll(i,t_plot))-min(UAll(i,t_plot));
  PP_U(i)=max(PAll(i,t_plot))-min(PAll(i,t_plot));
  %   RPSI_U(i)=max(diff(UAll(i,t_plot)))/mean(UAll(i,t_plot))*1e3;
  %   RPSI_P(i)=max(diff(PAll(i,t_plot)))/mean(PAll(i,t_plot))*1e3;
end

%% 计算锥度
% TaperRate=(Rin-Rout)/Rmean
TaperRate=(sqrt(MeanA(1)/pi)-sqrt(MeanA(3)/pi))/(sqrt(MeanA(1)/pi)+sqrt(MeanA(2)/pi)+sqrt(MeanA(3)/pi))*3;

%% Wave Intensity Analysis
figure;
PpAll=zeros(length(t_plot)-1,NumHisPt);
PnAll=PpAll;UpAll=PpAll;UnAll=PpAll;
for i=1:NumHisPt
  WIA_pt=i;
  VesID=1;
  [A beta]=Eval_beta_A(2*sqrt(MeanA(i)/1e12/pi),VesParam(4,VesID),VesParam(3,VesID));
  [Pp,Pn,Up,Un]=WIA(PAll(WIA_pt,t_plot,VesID)*133,UAll(WIA_pt,t_plot,VesID)/1e3,AAll(WIA_pt,t_plot,VesID)/1e12,beta,VesParam(4,:),VesParam(3,:),VesParam(6,:));
  Pp=Pp/133;Pn=Pn/133;Up=Up*1e3;Un=Un*1e3;
  PpAll(:,i)=Pp;
  PnAll(:,i)=Pn;
  UpAll(:,i)=Up;
  UnAll(:,i)=Un;
  % 前向波与反向波的PI
  PI_Pp(i)=(max(PpAll(:,i))-min(PpAll(:,i)))./mean(PpAll(:,i));
  PI_Pn(i)=(max(PnAll(:,i))-min(PnAll(:,i)))./mean(PnAll(:,i));
  PI_Up(i)=(max(UpAll(:,i))-min(UpAll(:,i)))./mean(UpAll(:,i));
  PI_Un(i)=(max(UnAll(:,i))-min(UnAll(:,i)))./mean(UnAll(:,i));
  % 归一化
  maxValue=max([PAll(WIA_pt,t_plot,VesID) max(PpAll(:,i)) max(PnAll(:,i))]);
  minValue=min([PAll(WIA_pt,t_plot,VesID) min(PpAll(:,i)) min(PnAll(:,i))]);
  NormP=(PAll(WIA_pt,t_plot,VesID)-minValue)./(maxValue-minValue);
  NormForwP=(PpAll(:,i)-minValue)./(maxValue-minValue);
  NormBackwP=(PnAll(:,i)-minValue)./(maxValue-minValue);
  subplot(1,NumHisPt,i);hold on;
  plot(tAll(2,t_plot),PAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
  plot(tAll(2,t_plot(1:end-1)),PpAll(:,i),'r','LineWidth',2);
  plot(tAll(2,t_plot(1:end-1)),PnAll(:,i),'k','LineWidth',2);
%   plot(tAll(2,t_plot),NormP,'b','LineWidth',2);
%   plot(tAll(2,t_plot(1:end-1)),NormForwP,'r--','LineWidth',2);
%   plot(tAll(2,t_plot(1:end-1)),NormBackwP,'k:','LineWidth',2);
%   axis([0 0.8 0 1]);
  xlabel('t(s)','FontName','Times New Roman','FontSize',20);
  ylabel('Norm. Pressure','FontName','Times New Roman','FontSize',20);
  set(gca,'FontName','Times New Roman','FontSize',20);
%   legend('Total','Forward','Backward','FontSize',20);
  switch i
    case 1
      title('始端','FontName','Times New Roman','FontSize',20);
    case 2
      title('中点','FontName','Times New Roman','FontSize',20);
    case 3
      title('末端','FontName','Times New Roman','FontSize',20);
  end
  
  % 归一化
  maxValue=max([UAll(WIA_pt,t_plot,VesID) max(UpAll(:,i)) max(UnAll(:,i))]);
  minValue=min([UAll(WIA_pt,t_plot,VesID) min(UpAll(:,i)) min(UnAll(:,i))]);
  NormU=(UAll(WIA_pt,t_plot,VesID)-minValue)./(maxValue-minValue);
  NormForwU=(UpAll(:,i)-minValue)./(maxValue-minValue);
  NormBackwU=(UnAll(:,i)-minValue)./(maxValue-minValue);
%   subplot(2,NumHisPt,NumHisPt+i);hold on;
  %   plot(UAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
  %   plot(UpAll(:,i),'r','LineWidth',2);
  %   plot(UnAll(:,i),'k','LineWidth',2);
%   plot(NormU,'b','LineWidth',2);
%   plot(NormForwU,'r','LineWidth',2);
%   %   NormBackwU=2*max(NormBackwU)-NormBackwU;
%   plot(NormBackwU,'k','LineWidth',2);
end

figure;
DampingRate=(PI_P/PI_P(1)-1)*100;
x=VesParam(1,:)/(VesParam(21,:)-1)*(0:VesParam(21,:)-1);
plot(x,DampingRate,'o-','LineWidth',2);
ylim([-20 20]);
xlabel('x(m)','FontName','Times New Roman','FontSize',20);
ylabel('PI/PI_{in}(%)','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);