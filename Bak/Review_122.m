%% Data review
close all;clc;
PeriodNum=4;    % Number of cardiac cycle in the simulation
DispPeriodNum=2;
colorOrder=64;
% Hemodynamic params calculation
t_plot=1000*Period*(PeriodNum-DispPeriodNum)+1:1000*Period*(PeriodNum);
all_plot=1000*Period*0+1:1000*Period*(PeriodNum);
[MeanP MeanQ MeanU MeanA MeanVisc tAll PAll UAll QAll AAll ViscAll WfAll WbAll...
  PforwAll PbackwAll UforwAll UbackwAll]=HemoAnal4Network(VesType,VesNum,t_plot,all_plot);
MeanP=MeanP/133;
MeanQ=MeanQ*1e12*60;
MeanU=MeanU*1000;
MeanA=MeanA*1e12;
MeanVisc=MeanVisc*1e3;
PAll=PAll/133;
QAll=QAll*1e12*60;
UAll=UAll*1000;
AAll=AAll*1e12;
ViscAll=ViscAll*1e3;
PforwAll=PforwAll/133;
PbackwAll=PbackwAll/133;
UforwAll=UforwAll*1000;
UbackwAll=UbackwAll*1000;

PI_U=zeros(VesNum,1);
PI_P=zeros(VesNum,1);
RI_U=zeros(VesNum,1);
RI_P=zeros(VesNum,1);
OrgPTT_U=zeros(VesNum,1);
OrgPTT_P=zeros(VesNum,1);
RPSI_U=zeros(VesNum,1);
RPSI_P=zeros(VesNum,1);

%% Calculate PI & PTT
for j=1:VesNum
  % PI
  PI_U(j)=(max(UAll(2,t_plot,j))-min(UAll(2,t_plot,j)))./MeanU(j,2);
  PI_P(j)=(max(PAll(3,t_plot,j))-min(PAll(3,t_plot,j)))./MeanP(j,3);
  % RI
  if mean(UAll(2,t_plot,j)>0)
    RI_U(j)=(max(UAll(2,t_plot,j))-min(UAll(2,t_plot,j)))./max(UAll(2,t_plot,j));
  else
    RI_U(j)=(max(UAll(2,t_plot,j))-min(UAll(2,t_plot,j)))./min(UAll(2,t_plot,j));
  end
  RI_P(j)=(max(PAll(2,t_plot,j))-min(PAll(2,t_plot,j)))./MeanP(j,2);
  % PTT of U and P
  % Foot-to-foot method
  if (mean(UAll(2,all_plot,j))>0 && j~=361)
    [signal,p_time,p_press,u_time,u_press]=findpu(UAll(3,all_plot,j)+10,1000,1);
    %         [signal,p_time,p_press,u_time,u_press]=findpu(UnAll(:,j),1000,0);
  else
    [signal,p_time,p_press,u_time,u_press]=findpu(UAll(3,all_plot,j),1000,0);
    %         [signal,p_time,p_press,u_time,u_press]=findpu(UnAll(:,j),1000,0);
  end
  indu=u_time(end);
  [signal,p_time,p_press,u_time,u_press]=findpu(PAll(3,all_plot,j),1000,1);
  %     [signal,p_time,p_press,u_time,u_press]=findpu(PnAll(:,j),1000,1);
  indp=u_time(end);
  % PTT for display
  %     OrgPTT_P(j)=round(indp+DispPeriodNum*Period*1000);
  %     OrgPTT_U(j)=round(indu+DispPeriodNum*Period*1000);
  OrgPTT_P(j)=round(indp);
  OrgPTT_U(j)=round(indu);
  % RPSI
  RPSI_U(j)=max(diff(UAll(2,t_plot,j)))/mean(UAll(2,t_plot,j))*1e3;
  RPSI_P(j)=max(diff(PAll(2,t_plot,j)))/mean(PAll(2,t_plot,j))*1e3;
end
% Process PI
PI_U(PI_U<0)=-PI_U(PI_U<0);
RI_U(RI_U<0)=-RI_U(RI_U<0);
PTT_P=OrgPTT_P-OrgPTT_P(1);
PTT_U=OrgPTT_U-OrgPTT_U(1);
% PTT_P=max(PTT_P)-PTT_P;
% PTT_U=max(PTT_U)-PTT_U;

% Produce the RGB components for plotting the distribution graph
PTT_U_RGB=ColorScal(PTT_U,colorOrder);
PTT_P_RGB=ColorScal(PTT_P,colorOrder);
PI_U_RGB=ColorScal(PI_U,colorOrder);
PI_P_RGB=ColorScal(PI_P,colorOrder);

%% Plot the distribution of PI & PTT
% PlotMorph(PI_U);title('PI_U');
% PlotMorph(PI_P);title('PI_P');
% % PlotMorph(RI_U);title('RI_U');
% % PlotMorph(RI_P);title('RI_P');
% PlotMorph(PTT_U);title('PTT_U');
% PlotMorph(PTT_P);title('PTT_P');

% Calculate order of the art & cap vessels
AllOrder=zeros(VesNum,1);
AllLenOrder=zeros(VesNum,1);
for i=1:VesNum
  tmpInd=i;
  order=1;
  len=Len(1);
  while 1
    if From(tmpInd)==1;
      AllOrder(i)=order;
      AllLenOrder(i)=len;
      break;
    else
      inInd=find(From(tmpInd)==To);
      if length(inInd)==1
        len=len+Len(tmpInd);
        tmpInd=inInd;
        order=order+1;
      elseif length(inInd)==2
        order=0;
        break;
      else
        order=0;
        break;
      end
    end
  end
end
AllOrder(AllOrder==0)=max(AllOrder)+1;
order_x=min(AllOrder):max(AllOrder);
order_x(end)=[];
order_x(1)=[];
% Plot the PTT & PI versus order
for i=1:length(order_x)
  ind=(AllOrder==order_x(i));
  PI_P_Plot(i)=mean(PI_P(ind));
  PI_P_Std(i)=std(PI_P(ind));
  PI_U_Plot(i)=mean(PI_U(ind));
  PI_U_Std(i)=std(PI_U(ind));
  %     RI_P_Plot(i)=mean(RI_P(ind));
  %     RI_U_Plot(i)=mean(RI_U(ind));
  PTT_P_Plot(i)=mean(PTT_P(ind));
  PTT_P_Std(i)=std(PTT_P(ind));
  PTT_U_Plot(i)=mean(PTT_U(ind));
  PTT_U_Std(i)=std(PTT_U(ind));
end
figure;hold on;whitebg('w');
errorbar(order_x,PI_P_Plot,PI_P_Std,'k.-','LineWidth',2,'MarkerSize',24);
errorbar(order_x,PI_U_Plot,PI_U_Std,'ko-','LineWidth',2,'MarkerSize',7);
xlabel('Generation number');
ylabel('PI of Pressure & Velocity waveforms');
legend('PI\_P','PI\_U');
% figure;hold on;whitebg('w');
% plot(order_x,RI_P_Plot,'k.-','LineWidth',2,'MarkerSize',24);
% plot(order_x,RI_U_Plot,'ko-','LineWidth',2,'MarkerSize',7);
% xlabel('Generation number');
% ylabel('RI of Pressure & Velocity waveforms');
% legend('RI\_P','RI\_U');
figure;hold on;whitebg('w');
errorbar(order_x,PTT_P_Plot,PTT_P_Std,'k.-','LineWidth',2,'MarkerSize',24);
errorbar(order_x,PTT_U_Plot,PTT_U_Std,'ko-','LineWidth',2,'MarkerSize',7);
xlabel('Generation number');
ylabel('PTT of Pressure & Velocity waveforms (ms)');
legend('PTT\_P','PTT\_U');

%% Plot the PI & PTT versus len
len_x=linspace(AllLenOrder(1),max(AllLenOrder),10);
for i=1:length(len_x)-1
  ind=(AllLenOrder<len_x(i+1) & AllLenOrder>=len_x(i));
  PI_P_Plot(i)=mean(PI_P(ind));
  PI_P_Std(i)=std(PI_P(ind));
  PI_U_Plot(i)=mean(PI_U(ind));
  PI_U_Std(i)=std(PI_U(ind));
  PTT_P_Plot(i)=mean(PTT_P(ind));
  PTT_P_Std(i)=std(PTT_P(ind));
  PTT_U_Plot(i)=mean(PTT_U(ind));
  PTT_U_Std(i)=std(PTT_U(ind));
end
figure;hold on;whitebg('w');
subplot(121)
errorbar(PI_P_Plot(1:length(len_x)-1),PI_P_Std(1:length(len_x)-1),'k.-','LineWidth',2,'MarkerSize',24);
% errorbar(order_x,PI_U_Plot,PI_U_Std,'ko-','LineWidth',2,'MarkerSize',7);
xlabel('Accumulated Length');
ylabel('PI_P');
legend('PI\_P');
% figure;hold on;whitebg('w');
subplot(122);hold on;
errorbar(PTT_P_Plot(1:length(len_x)-1),PTT_P_Std(1:length(len_x)-1),'k.-','LineWidth',2,'MarkerSize',24);
% errorbar(order_x,PTT_U_Plot,PTT_U_Std,'ko-','LineWidth',2,'MarkerSize',7);
xlabel('Accumulated Length');
ylabel('PTT_P(ms)');
legend('PTT\_P');

%% Plot the PTT & PI versus Vessel Type
for i=1:3
  tmpInd=find(VesCategory==i);
  MeanPI_P(i)=mean(PI_P(tmpInd));
  StdPI_P(i)=std(PI_P(tmpInd));
  CVPI_P(i)=StdPI_P(i)/MeanPI_P(i);
  MeanPI_U(i)=mean(PI_U(tmpInd));
  StdPI_U(i)=std(PI_U(tmpInd));
  CVPI_U(i)=StdPI_U(i)/MeanPI_U(i);
  MeanPTT_P(i)=mean(PTT_P(tmpInd));
  StdPTT_P(i)=std(PTT_P(tmpInd));
  CVPTT_P(i)=StdPTT_P(i)/MeanPTT_P(i);
  MeanPTT_U(i)=mean(PTT_U(tmpInd));
  StdPTT_U(i)=std(PTT_U(tmpInd));
  CVPTT_U(i)=StdPTT_U(i)/MeanPTT_U(i);
  MeanRPSI_P(i)=mean(RPSI_P(tmpInd));
  StdRPSI_P(i)=std(RPSI_P(tmpInd));
  CVRPSI_P(i)=StdRPSI_P(i)/MeanRPSI_P(i);
end
tmpVesType=VesCategory;
tmpVesType(MeanPTTInd)=[];
PTT_U(MeanPTTInd)=[];
for i=1:3
    tmpInd=find(tmpVesType==i);
    MeanPTT_U(i)=mean(PTT_U(tmpInd));
    StdPTT_U(i)=std(PTT_U(tmpInd));
    CVPTT_U(i)=StdPTT_U(i)/MeanPTT_U(i);
end
figure;hold on;whitebg('w');
errorbar(MeanPI_P,StdPI_P,'k.-','MarkerSize',24);title('PI_P vs. Vessel Type');
figure;hold on;whitebg('w');
errorbar(MeanPI_U,StdPI_U,'k.-','MarkerSize',24);title('PI_U vs. Vessel Type');
figure;hold on;whitebg('w');
errorbar(MeanPTT_P,StdPTT_P,'k.-','MarkerSize',24);title('PTT_P vs. Vessel Type');
figure;hold on;whitebg('w');
errorbar(MeanPTT_U,StdPTT_U,'k.-','MarkerSize',24);title('PTT_U vs. Vessel Type');

AllOrder(AllOrder==0)=max(AllOrder)+1;
AllOrder=max(AllOrder)-AllOrder+1;
Order_RGB=ColorScal(AllOrder,256);

%% Plot the hemodynamic parameters of each vessel
for j=1:VesNum
  h=figure;
  subplot(221);hold on;
  %     plot(tAll(1,t_plot,j), PAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll(2,t_plot,j), 'b', 'LineWidth', 2);
  plot(tAll(2,round(OrgPTT_P(j)),j), PAll(2,round(OrgPTT_P(j)),j), 'ro');
  %     plot(tAll(2,all_plot(t_plot(2):end),j), PpAll(t_plot(1):end,j), 'r', 'LineWidth', 2);
  %     plot(tAll(2,all_plot(t_plot(2):end),j), PnAll(t_plot(1):end,j), 'k', 'LineWidth', 2);
  %     plot(tAll(3,t_plot,j), PAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Pressure(mmHg)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Pressure Vessel ' num2str(SegName(j),'%04d')]);
  
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
  plot(tAll(2,round(OrgPTT_U(j)),j), UAll(2,round(OrgPTT_U(j)),j), 'ro');
  %     plot(tAll(2,all_plot(t_plot(2):end),j), UpAll(t_plot(1):end,j), 'r', 'LineWidth', 2);
  %     plot(tAll(2,all_plot(t_plot(2):end),j), UnAll(t_plot(1):end,j), 'k', 'LineWidth', 2);
  %     plot(tAll(3,t_plot,j), UAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Velocity(mm/s)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Velocity Vessel ' num2str(SegName(j),'%04d')]);
  
  subplot(223);hold on;
  plot(tAll(2,t_plot,j), ViscAll(2,t_plot,j), 'b', 'LineWidth', 2);
  ylabel('Viscosity(mPam)');xlabel('t(s)');
  title(['Viscosity Vessel ' num2str(SegName(j),'%04d')]);
  
  subplot(224);hold on;
  plot(tAll(2,t_plot,j), PAll(2,t_plot,j), 'b', 'LineWidth', 2);
  ylabel('Area(um2)');xlabel('t(s)');
  title(['Area Vessel ' num2str(SegName(j),'%04d')]);
  
  %     saveas(h,['Vessel ' num2str(SegName(j),'%04d') '_' num2str(j,'%03d')]);
  %     saveas(h,['Vessel ' num2str(SegName(j),'%04d') '_' num2str(j,'%03d')],'png');
  close(h);
end

%% WIA
for j=1:VesNum
  h=figure;
  for i=1:3
    PpAll=zeros(length(t_plot)-1,1);
    PnAll=PpAll;UpAll=PpAll;UnAll=PpAll;
    WIA_pt=i;
    VesID=j;
    [A beta]=Eval_beta_A(VesParam(2,VesID),VesParam(4,VesID),VesParam(3,VesID));
    [Pp,Pn,Up,Un]=WIA(PAll(WIA_pt,t_plot,VesID)*133,UAll(WIA_pt,t_plot,VesID)/1e3,AAll(WIA_pt,t_plot,VesID)/1e12,beta,VesParam(4,:),VesParam(3,:),VesParam(6,:));
    Pp=Pp/133;Pn=Pn/133;Up=Up*1e3;Un=Un*1e3;
    PpAll=Pp;
    PnAll=Pn;
    UpAll=Up;
    UnAll=Un;
    subplot(2,3,i);hold on;
    plot(PAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
    plot(PpAll,'r','LineWidth',2);
    plot(PnAll,'k','LineWidth',2);
    subplot(2,3,i+3);hold on;
    plot(UAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
    plot(UpAll,'r','LineWidth',2);
    plot(UnAll,'k','LineWidth',2);
  end
  saveas(h,['VesselWIA ' num2str(SegName(j),'%04d') '_' num2str(j,'%03d')],'png');
  close(h);
end

%% Compare the 1D simulation data with data provided
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
% plot(FuncPara(:,7),'r.-');
if MeasDiam
  plot(DataArray(:,8),'r.-');
else
  plot(FuncPara(:,7),'r.-');
end
MeanU(MeanU(:,2)<0,2)=-MeanU(MeanU(:,2)<0,2);
plot(MeanU(:,2),'.-');
ylabel('Flow velocity (mm/s)');
legend('Pries data','1D data');
subplot(224);hold on;
if MeasDiam
  plot(DataArray(:,5).*DataArray(:,5)*pi/4,'r.-');
else
  plot(FuncPara(:,3).*FuncPara(:,3)*pi/4,'r.-');
end
plot(MeanA(:,2),'.-');
ylabel('Area (um2)');
legend('Pries data','1D data');

%% Compare waveforms of different vessels
figure;hold on;
t_plot=1:1000*Period*PeriodNum;
Ind=[1 2001];
tmpInd1=find(SegName==Ind(1));
tmpInd2=find(SegName==Ind(2));
% disp(SegName(Ind(2)));
if max(tmpInd1)>546
  fprintf(stderr,'Maximum vessel number is 546\n');
else
  % Scale first
  %     Tmp1=(UAll(1,t_plot,tmpInd1)-mean(UAll(1,t_plot,tmpInd1)))/(max(UAll(1,t_plot,tmpInd1))-min(UAll(1,t_plot,tmpInd1)));
  %     Tmp2=(UAll(1,t_plot,tmpInd2)-mean(UAll(1,t_plot,tmpInd2)))/(max(UAll(1,t_plot,tmpInd2))-min(UAll(1,t_plot,tmpInd2)));
  %     Tmp1=UAll(1,t_plot,tmpInd1)/max(UAll(1,t_plot,tmpInd1));
  %     Tmp2=UAll(1,t_plot,tmpInd2)/max(UAll(1,t_plot,tmpInd2));
  Tmp1=UAll(1,t_plot,tmpInd1);
  Tmp2=UAll(1,t_plot,tmpInd2);
  % Plot scaled curves
  plot(tAll(1,t_plot,tmpInd1),Tmp1,'.-');
  plot(tAll(1,t_plot,tmpInd2),Tmp2,'r.-');
end
