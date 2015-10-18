%% Data review of the symmetric network
close all;clc;
% Hemodynamic params calculation
PeriodNum=6;    % Number of cardiac cycle in the simulation
% t_plot=1:800*PeriodNum;
t_plot=800*(PeriodNum-1)+1:800*PeriodNum;
all_plot=1:800*PeriodNum;
NumHisPt=VesParam(21,1);
[MeanP MeanQ MeanU MeanA MeanVisc tAll PAll UAll QAll AAll ViscAll]...
  =HemoAnal4Network(NumHisPt,VesType,VesNum,t_plot,all_plot);

% Tree order=5
% R1=8*MeanVisc(1,2).*LenVec(1)/pi./(((MeanA(1,2)+MeanA(1,2))/2/pi).^2);
% R2=8*MeanVisc(2,2).*LenVec(2)/pi./(((MeanA(2,2)+MeanA(2,2))/2/pi).^2);
% R3=8*MeanVisc(4,2).*LenVec(4)/pi./(((MeanA(4,2)+MeanA(3,2))/2/pi).^2);
% R4=8*MeanVisc(8,2).*LenVec(8)/pi./(((MeanA(8,2)+MeanA(8,2))/2/pi).^2);
% R5=8*MeanVisc(16,2).*LenVec(16)/pi./(((MeanA(16,2)+MeanA(16,2))/2/pi).^2);
% R=R1+R2/2+R3/4+R4/8+R5/16;
% dA1=(MeanA(1,1)-MeanA(1,3))/mean(MeanA(1,:));
% dA2=(MeanA(2,1)-MeanA(2,3))/mean(MeanA(2,:));
% dA3=(MeanA(4,1)-MeanA(4,3))/mean(MeanA(4,:));
% dA4=(MeanA(8,1)-MeanA(8,3))/mean(MeanA(8,:));
% dA5=(MeanA(16,1)-MeanA(16,3))/mean(MeanA(16,:));
% dAAll=[dA1 dA2 dA3 dA4 dA5];
% dA=mean(dAAll);

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
  
%   print(h,'-dpng',['Vessel ' num2str(SegName(j),'%04d')]);
  close(h);
end

% Calculate PI
PI_U=zeros(VesNum,3);
PI_P=zeros(VesNum,3);
RPSI_U=zeros(VesNum,3);
RPSI_P=zeros(VesNum,3);
for i=1:VesNum
  for j=1:VesParam(21,1) % start, middle, end point
    PI_U(i,j)=(max(UAll(j,t_plot,i))-min(UAll(j,t_plot,i)))/mean(UAll(j,t_plot,i));
    PI_P(i,j)=(max(PAll(j,t_plot,i))-min(PAll(j,t_plot,i)))/mean(PAll(j,t_plot,i));
    % RPSI
    RPSI_U(i,j)=max(diff(UAll(j,t_plot,i)))/mean(UAll(j,t_plot,i))*1e3;
    RPSI_P(i,j)=max(diff(PAll(j,t_plot,i)))/mean(PAll(j,t_plot,i))*1e3;
  end
end

%% WIA
figure;
PpAll=zeros(3,length(t_plot)-1,3);
PnAll=PpAll;UpAll=PpAll;UnAll=PpAll;
cnt=0;
for j=1:1
% for j=[1 2 4 8 16]
  cnt=cnt+1;
  for i=1:3
    WIA_pt=i;
    VesID=j;
    [A beta]=Eval_beta_A(2*sqrt(MeanA(j,i)/1e12/pi),VesParam(4,VesID),VesParam(3,VesID));
    [Pp,Pn,Up,Un]=WIA(PAll(WIA_pt,t_plot,VesID)*133,UAll(WIA_pt,t_plot,VesID)/1e3,AAll(WIA_pt,t_plot,VesID)/1e12,beta,VesParam(4,:),VesParam(3,:),VesParam(6,:));
    Pp=Pp/133;Pn=Pn/133;Up=Up*1e3;Un=Un*1e3;
    PpAll(WIA_pt,:,VesID)=Pp;
    PnAll(WIA_pt,:,VesID)=Pn;
    UpAll(WIA_pt,:,VesID)=Up;
    UnAll(WIA_pt,:,VesID)=Un;
    maxValue=max([PAll(WIA_pt,t_plot,VesID) max(PpAll(WIA_pt,:,VesID)) max(PnAll(WIA_pt,:,VesID))]);
    minValue=min([PAll(WIA_pt,t_plot,VesID) min(PpAll(WIA_pt,:,VesID)) min(PnAll(WIA_pt,:,VesID))]);
    NormP=(PAll(WIA_pt,t_plot,VesID)-minValue)./(maxValue-minValue);
    NormForwP=(PpAll(WIA_pt,:,VesID)-minValue)./(maxValue-minValue);
    NormBackwP=(PnAll(WIA_pt,:,VesID)-minValue)./(maxValue-minValue);
%     subplot(3,VesNum,cnt+(i-1)*VesNum);hold on;
    subplot(1,3,i);hold on;
%     plot(tAll(1,t_plot,1),PAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
%     plot(tAll(1,t_plot(1:end-1),1),PpAll(WIA_pt,:,VesID),'r--','LineWidth',2);
%     plot(tAll(1,t_plot(1:end-1),1),PnAll(WIA_pt,:,VesID),'k:','LineWidth',2);
    plot(tAll(2,t_plot,2),NormP,'b','LineWidth',2);
    plot(tAll(2,t_plot(1:end-1),2),NormForwP,'r--','LineWidth',2);
    plot(tAll(2,t_plot(1:end-1),2),NormBackwP,'k:','LineWidth',2);
    xlabel('t(s)','FontName','Times New Roman','FontSize',20);
    ylabel('Norm. Pressure','FontName','Times New Roman','FontSize',20);
    set(gca,'FontName','Times New Roman','FontSize',20);
    legend('Total','Forward','Backward','FontSize',20);
    switch i
      case 1
        title('Ê¼¶Ë','FontName','Times New Roman','FontSize',20);
      case 2
        title('ÖÐµã','FontName','Times New Roman','FontSize',20);
      case 3
        title('Ä©¶Ë','FontName','Times New Roman','FontSize',20);
    end
%     subplot(3,VesNum,j+(i+2)*VesNum);hold on;
%     plot(UAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
%     plot(UpAll(WIA_pt,:,VesID),'r','LineWidth',2);
%     plot(UnAll(WIA_pt,:,VesID),'k','LineWidth',2);
  end
end

%% Reflection coefficient
rho=1050;   % kg/m3
[A beta]=Eval_beta_A(VesParam(2,:),VesParam(4,:),VesParam(3,:));
c=sqrt(beta./2./rho).*MeanA(:,2)'.^0.25;
Rf=(MeanA(1,2)/c(1)-MeanA(2,2)/c(2)-MeanA(3,2)/c(3))/(MeanA(1,2)/c(1)+MeanA(2,2)/c(2)+MeanA(3,2)/c(3));

(PI_P(1,1)-PI_P(end,3))/PI_P(1,1)
(PI_U(1,1)-PI_U(end,3))/PI_U(1,1)

%% Damping Rate
figure;
DampingRate=(PI_P(1,:)/PI_P(1,1)-1)*100;
x=VesParam(1,1)/(VesParam(21,1)-1)*(0:VesParam(21,1)-1);
plot(x,DampingRate,'o-','LineWidth',2);
ylim([-30 30]);
xlabel('x(m)','FontName','Times New Roman','FontSize',20);
ylabel('PI/PI_{in}(%)','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);
