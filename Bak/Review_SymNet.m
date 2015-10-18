%% Data review of the symmetric network
close all;clc;
% Hemodynamic params calculation
PeriodNum=6;    % Number of cardiac cycle in the simulation
t_plot=800*(PeriodNum-1)+1:800*PeriodNum;
all_plot=1:800*PeriodNum;
NumHisPt=VesParam(21,1);
[MeanP MeanQ MeanU MeanA MeanVisc tAll PAll UAll QAll AAll ViscAll]...
  =HemoAnal4Network(NumHisPt,VesType,VesNum,t_plot,all_plot);
R1=8*MeanVisc(1,2).*LenVec(1)/pi./(((MeanA(1,2)+MeanA(46,2))/2/pi).^2);
R2=8*MeanVisc(2,2).*LenVec(2)/pi./(((MeanA(2,2)+MeanA(45,2))/2/pi).^2);
R3=8*MeanVisc(4,2).*LenVec(4)/pi./(((MeanA(4,2)+MeanA(43,2))/2/pi).^2);
R4=8*MeanVisc(8,2).*LenVec(8)/pi./(((MeanA(8,2)+MeanA(39,2))/2/pi).^2);
R5=8*MeanVisc(16,2).*LenVec(16)/pi./(((MeanA(16,2)+MeanA(31,2))/2/pi).^2);
R=R1*2+R2+R3/2+R4/4+R5/16;
RAll=(MeanP(1,1)-MeanP(46,3))/MeanQ(46,3);
dA=(MeanA(1,1)-MeanA(46,3))/MeanA(1,1);

MeanP=MeanP/133;
MeanQ=MeanQ*1e12*60;
MeanU=MeanU*1000;
MeanA=MeanA*1e12;
PAll=PAll/133;
QAll=QAll*1e12*60;
UAll=UAll*1000;
AAll=AAll*1e12;

% Each vessel's P,Q,U
for j=1:VesNum
  h=figure;
  subplot(221);hold on;
  plot(tAll(1,t_plot,j), PAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll(2,t_plot,j), 'b', 'LineWidth', 2);
  plot(tAll(3,t_plot,j), PAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Pressure(mmHg)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Pressure Vessel ' num2str(SegName(j),'%04d')]);
  
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
  
  subplot(224);hold on;
  plot(tAll(1,t_plot,j), PAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll(2,t_plot,j), 'b', 'LineWidth', 2);
  plot(tAll(3,t_plot,j), PAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Area(um2)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Area Vessel ' num2str(SegName(j),'%04d')]);
  
  %     print(h,'-dpng',['Vessel ' num2str(SegName(j),'%04d')]);
  close(h);
end

% Calculate PI
PI_U=zeros(46,3);
PI_P=zeros(46,3);
RPSI_U=zeros(46,3);
RPSI_P=zeros(46,3);
PP=zeros(46,3);
P_Mean=zeros(46,3);
for i=1:VesNum
  for j=1:3
    % PI
    PI_U(i,j)=(max(UAll(j,t_plot,i))-min(UAll(j,t_plot,i)))/mean(UAll(j,t_plot,i));
    PI_P(i,j)=(max(PAll(j,t_plot,i))-min(PAll(j,t_plot,i)))/mean(PAll(j,t_plot,i));
    % RPSI
    RPSI_U(i,j)=max(diff(UAll(j,t_plot,i)))/mean(UAll(j,t_plot,i))*1e3;
    RPSI_P(i,j)=max(diff(PAll(j,t_plot,i)))/mean(PAll(j,t_plot,i))*1e3;
%     % PP
%     PP(i,j)=(max(PAll(j,t_plot,i))-min(PAll(j,t_plot,i)));
%     % P_Mean
%     P_Mean(i,j)=mean(PAll(j,t_plot,i));
  end
end

PI_Reduce_Ratio=(PI_P(1,1)-PI_P(46,3))/PI_P(1,1)
RPSI_Reduce_Ratio=(RPSI_P(1,1)-RPSI_P(46,3))/RPSI_P(1,1)

%% WIA
PI_Pp=zeros(9,3);
PN_Ratio=zeros(9,3);
WIA_Ind=[1 2 4 8 16 32 40 44 46];
figure;
for j=1:length(WIA_Ind)
  for i=1:3
    PpAll=zeros(length(t_plot)-1,1);
    PnAll=PpAll;UpAll=PpAll;UnAll=PpAll;
    WIA_pt=i;
    VesID=WIA_Ind(j);
    [A beta]=Eval_beta_A(VesParam(2,VesID),VesParam(4,VesID),VesParam(3,VesID));
    [Pp,Pn,Up,Un]=WIA(PAll(WIA_pt,t_plot,VesID)*133,UAll(WIA_pt,t_plot,VesID)/1e3,AAll(WIA_pt,t_plot,VesID)/1e12,beta,VesParam(4,:),VesParam(3,:),VesParam(6,:));
    Pp=Pp/133;Pn=Pn/133;Up=Up*1e3;Un=Un*1e3;
    PpAll=Pp;
    PnAll=Pn;
    UpAll=Up;
    UnAll=Un;
    PI_Pp(j,i)=(max(PpAll)-min(PpAll))/mean(PpAll);
    PN_Ratio(j,i)=max(PpAll)/max(PnAll);
    % PP
    PP(j,i)=(max(PAll(WIA_pt,t_plot,VesID))-min(PAll(WIA_pt,t_plot,VesID)));
    % P_Mean
    P_Mean(j,i)=mean(PAll(WIA_pt,t_plot,VesID));
%     figure;hold on
    subplot(6,length(WIA_Ind),j+(i-1)*length(WIA_Ind));hold on;
    plot(PAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
    plot(PpAll,'r','LineWidth',2);
    plot(PnAll,'k','LineWidth',2);
    subplot(6,length(WIA_Ind),j+(i+2)*length(WIA_Ind));hold on;
    plot(UAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
    plot(UpAll,'r','LineWidth',2);
    plot(UnAll,'k','LineWidth',2);
  end
end

% Compare the 1D simulation data with data provided
figure;
subplot(221);hold on;
plot(MeanP(:,2),'.-');
ylabel('Pressure(mmHg)');
legend('1D data');
subplot(222);hold on;
plot(MeanQ(:,2),'.-');
ylabel('Flow rate(nl/min)');
legend('1D data');
subplot(223);hold on;
plot(MeanU(:,2),'.-');
ylabel('Flow velocity (mm/s)');
legend('1D data');
subplot(224);hold on;
plot(MeanA(:,2),'.-');
ylabel('Area (um2)');
legend('1D data');

% Compare waveforms of different vessels
figure;hold on;
% t_plot=1:800*PeriodNum;
Ind=[1 22];
tmpInd1=find(SegName==Ind(1));
tmpInd2=find(SegName==Ind(2));
% disp(SegName(Ind(2)));
if max(tmpInd1)>546
  fprintf(stderr,'Maximum vessel number is 546\n');
else
  % Scale first
  Tmp1=(UAll(1,t_plot,tmpInd1)-mean(UAll(1,t_plot,tmpInd1)))/(max(UAll(1,t_plot,tmpInd1))-min(UAll(1,t_plot,tmpInd1)));
  Tmp2=(UAll(1,t_plot,tmpInd2)-mean(UAll(1,t_plot,tmpInd2)))/(max(UAll(1,t_plot,tmpInd2))-min(UAll(1,t_plot,tmpInd2)));
  % Plot scaled curves
  plot(tAll(1,t_plot,tmpInd1),Tmp1,'.-');
  plot(tAll(1,t_plot,tmpInd2),Tmp2,'r.-');
end

plot_pt=2;
figure;
% For order 6
% subplot(121)
% hold on
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,1),'b.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,2),'r.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,4),'g.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,8),'k.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,16),'m.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,32),'c.-')
% subplot(122)
% hold on
% ylim([0 1.8]);
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,94),'b.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,92),'r.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,88),'g.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,80),'k.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,64),'m.-')
% plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,32),'c.-')
% AllLen=LenVec(1)+LenVec(2)+LenVec(4)+LenVec(8)+LenVec(16)+LenVec(32);

% For order 5
subplot(221)
hold on
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,1),'b.-')
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,2),'r.-')
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,4),'g.-')
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,8),'k.-')
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,16),'m.-')
xlabel('t(s)');
ylabel('Velocity(mm/s)');
subplot(222)
hold on
% ylim([0 1.8]);
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,46),'b.-')
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,44),'r.-')
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,40),'g.-')
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,32),'k.-')
plot(tAll(1,t_plot,1),UAll(plot_pt,t_plot,16),'m.-')
xlabel('t(s)');
ylabel('Velocity(mm/s)');
subplot(223)
hold on
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,1),'b.-')
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,2),'r.-')
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,4),'g.-')
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,8),'k.-')
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,16),'m.-')
xlabel('t(s)');
ylabel('Pressure(mmHg)');
subplot(224)
hold on
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,46),'b.-')
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,44),'r.-')
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,40),'g.-')
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,32),'k.-')
plot(tAll(1,t_plot,1),PAll(plot_pt,t_plot,16),'m.-')
xlabel('t(s)');
ylabel('Pressure(mmHg)');
AllLen=LenVec(1)+LenVec(2)+LenVec(4)+LenVec(8)+LenVec(16);
% AllR=R1(1)+R1(2)+R1(4)+R1(8)+R1(16)+R1(32)+R1(40)+R1(44)+R1(46);