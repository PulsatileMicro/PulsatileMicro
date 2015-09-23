%% Egg818 Network Results Review
% Before running this script, what you have to do are:
% 1. Run main_egg.m (in the directory with main_546.m file)
% 2. Enter the data directory (such as 546_qR_E1)
% 3. Run this script and click "add to path" but not "change directory"
close all;clc;

%% Load data
saveFlag=1; % The flag of saving the workspace
list = ls;  % List the files in the directory
[Row Col]=size(list);
for i=1:Row
  % If Egg818.mat exist, load the data from .mat file but not .his files
  if strfind(list(i,:),'Egg818.mat')
    load Egg818.mat;
    saveFlag=0;
    break;
  end
end
% If .mat file doesn't exist, load data from .his files
if saveFlag
  PeriodNum=6;    % Number of cardiac cycles in the simulation
  DispPeriodNum=2;  % Number of cycles for display
  colorOrder=64;
  % Hemodynamic params calculation
  t_plot=1000*Period*(PeriodNum-DispPeriodNum)+1:1000*Period*(PeriodNum);
  all_plot=1000*Period*0+1:1000*Period*(PeriodNum);
  NumHisPt=VesParam(21,1);
  [MeanP MeanQ MeanU MeanA MeanVisc tAll PAll UAll QAll AAll ViscAll]...
    =HemoAnal4Network(NumHisPt,VesType,VesNum,t_plot,all_plot);
  MeanP=MeanP/133;
  MeanQ=MeanQ*1e12*60;
  MeanU=MeanU*1000;
  MeanA=MeanA*1e12;
  PAll=PAll/133;
  QAll=QAll*1e12*60;
  UAll=UAll*1000;
  AAll=AAll*1e12;
  PI_U=zeros(VesNum,1);
  PI_P=zeros(VesNum,1);
  OrgPTT_U=zeros(VesNum,1);
  OrgPTT_P=zeros(VesNum,1);
  RPSI_U=zeros(VesNum,1);
  RPSI_P=zeros(VesNum,1);
  SDRatio_P=zeros(VesNum,1);
  SDRatio_U=zeros(VesNum,1);
  saveFlag=1;
end

%% Calculate PI, PTT, RPSI
% PI & RPSI is computed based on the waveform in the middle point of each segment
% PTT is computed at the beginning of each segment
for j=1:VesNum
  % PI
  PI_U(j)=(max(UAll(2,t_plot,j))-min(UAll(2,t_plot,j)))./MeanU(j,2);
  PI_P(j)=(max(PAll(2,t_plot,j))-min(PAll(2,t_plot,j)))./MeanP(j,2);
  % PI_P(j)=(max(PAll(2,t_plot,j))-min(PAll(2,t_plot,j)))./(min(PAll(2,t_plot,j))+1/3*(max(PAll(2,t_plot,j))-min(PAll(2,t_plot,j))));
  
  % SDRatio
  SDRatio_U(j)=max(UAll(2,t_plot,j))/min(UAll(2,t_plot,j));
  % PTT of U and P
  % Foot-to-foot method
  if (mean(UAll(2,all_plot,j))>0 && j~=361)
%     [signal,p_time,p_press,u_time,u_press]=findpu(UAll(1,all_plot,j)+10,1000,1);
  else
%     [signal,p_time,p_press,u_time,u_press]=findpu(UAll(1,all_plot,j),1000,0);
  end
%   indu=u_time(end);
%   [signal,p_time,p_press,u_time,u_press]=findpu(PAll(1,all_plot,j),1000,1);
%   indp=u_time(end);
  
  % PTT for display
%   OrgPTT_P(j)=round(indp);
%   OrgPTT_U(j)=round(indu);
  
  % RPSI
  tmpU1=UAll(2,t_plot(1:20:end),j);
  tmpU=(tmpU1(1:end-2)+tmpU1(2:end-1)+tmpU1(3:end))/3;
  RPSI_U(j)=max(diff(tmpU)/20)/mean(tmpU)*1000;
%   RPSI_U(j)=max(diff(UAll(2,t_plot,j)))/mean(UAll(2,t_plot,j))*1e3;
%   RPSI_P(j)=max(diff(PAll(2,t_plot,j)))/mean(PAll(2,t_plot,j))*1e3;
end

ArtInd=[4 11 136 145 228 79 81 107 94 158 168 178 180 186 1026 922 905 740 736 794 791 820 669];
VenInd=[434 351 338 312 577 573 567 537 470 462 500 393 381 382 298 301 393 689 309 531 542 502 505];

ArtRPSI=[];VenRPSI=[];
for i=1:VesNum
  if ~isempty(find(SegName(i)==ArtInd))
%     SegName(i)
    ArtRPSI=[ArtRPSI [RPSI_U(i);SegName(i);MeanU(i,2);Diam(i)]];
  elseif ~isempty(find(SegName(i)==VenInd))
    SegName(i)
    VenRPSI=[VenRPSI [RPSI_U(i);SegName(i);MeanU(i,2);Diam(i)]];
  end
end

% Generate the RGB components for plotting the distribution map
PTT_P=OrgPTT_P-OrgPTT_P(1);
PTT_U=OrgPTT_U-OrgPTT_U(1);
PTT_U_RGB=ColorScal(PTT_U,colorOrder);
PTT_P_RGB=ColorScal(PTT_P,colorOrder);
PI_U_RGB=ColorScal(PI_U,colorOrder);
PI_P_RGB=ColorScal(PI_P,colorOrder);
RPSI_U_RGB=ColorScal(RPSI_U,colorOrder);
RPSI_P_RGB=ColorScal(RPSI_P,colorOrder);

%% Plot the distribution of PI,PTT,RPSI
% PlotMorph_Egg818(PI_U);title('PI_U');
PlotMorph_Egg818(PI_P);title('PI_P');
% PlotMorph_Egg818(PTT_U);title('PTT_U');
% PlotMorph_Egg818(PTT_P);title('PTT_P');
% PlotMorph_Egg818(RPSI_U);title('RPSI_U');
% PlotMorph_Egg818(RPSI_P);title('RPSI_P');

%% Compute order and accumulated length of the art & cap vessels
AllOrder=zeros(VesNum,1);
AllLenOrder=zeros(VesNum,1);
for i=1:VesNum
  tmpInd=i;
  order=1;
  len=Len(1);
  while 1
    if From(tmpInd)==830;
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

%% Plot the PI & PTT vs. Length
len_x=linspace(AllLenOrder(1),max(AllLenOrder),10); % 分为10等份
for i=1:length(len_x)-1
  ind=(AllLenOrder<len_x(i+1) & AllLenOrder>=len_x(i));
  %     i
  %     find(ind)
  MeanPIP_ord(i)=mean(PI_P(ind));
  StdPIP_ord(i)=std(PI_P(ind));
  MeanPIU_ord(i)=mean(PI_U(ind));
  StdPIU_ord(i)=std(PI_U(ind));
  MeanPTTP_ord(i)=mean(PTT_P(ind));
  StdPTTP_ord(i)=std(PTT_P(ind));
  MeanPTTU_ord(i)=mean(PTT_U(ind));
  StdPTTU_ord(i)=std(PTT_U(ind));
end
figure;hold on;whitebg('w');
subplot(121)
errorbar(MeanPIP_ord(1:length(len_x)-1),StdPIP_ord(1:length(len_x)-1),'k.-','LineWidth',2,'MarkerSize',24);
% errorbar(order_x,MeanPIU_ord,StdPIU_ord,'ko-','LineWidth',2,'MarkerSize',7);
xlabel('Accumulated Length');
ylabel('PI_P');
legend('PI\_P');
subplot(122);
errorbar(MeanPTTP_ord(1:length(len_x)-1),StdPTTP_ord(1:length(len_x)-1),'k.-','LineWidth',2,'MarkerSize',24);
% errorbar(order_x,MeanPTTU_ord,StdPTTU_ord,'ko-','LineWidth',2,'MarkerSize',7);
xlabel('Accumulated Length');
ylabel('PTT_P(ms)');
legend('PTT\_P');

%% PP vs. Pressure, PI vs. Pressure
figure;
subplot(121)
[r m b]=regression(MeanP(:,2)',PI_P'.*MeanP(:,2)');
p_x=0:0.01:90;
p_y=p_x.*m+b;
plot(MeanP(:,2),PI_P.*MeanP(:,2),'k.');
hold on;
plot(p_x,p_y,'LineWidth',2.5);
axis([0 90 0 90]);
xlabel('Mean Pressure(mmHg)');
ylabel('Pressure Amplitude(mmHg)');
subplot(122)
plot(MeanP(:,2),PI_P,'ko');
xlabel('Mean Pressure(mmHg)');
ylabel('PI');

%% Plot the PTT & PI vs. Vessel Type
R=128*Visc.*Len/pi./(Diam.^4);
for i=1:3
  tmpInd=find(VesCategory==i);
%   if i==1
%      tmpInd(SegName(tmpInd)>200 & SegName(tmpInd)<400)=[];
%   end
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
  MeanSDRatio_U(i)=mean(SDRatio_U(tmpInd));
  StdSDRatio_U(i)=std(SDRatio_U(tmpInd));
  
  % Study the morphological difference between types
  Diam_Typ(i)=mean(Diam(tmpInd));
  Diam_max(i)=max(Diam(tmpInd));
  Diam_min(i)=min(Diam(tmpInd));
  R_Typ(i)=mean(R(tmpInd));
end
figure;hold on;whitebg('w');
subplot(221)
errorbar(MeanPI_P,StdPI_P,'k.-','MarkerSize',24);title('PI_P vs. Vessel Type');
subplot(222)
errorbar(MeanPI_U,StdPI_U,'k.-','MarkerSize',24);title('PI_U vs. Vessel Type');
subplot(223)
errorbar(MeanPTT_P,StdPTT_P,'k.-','MarkerSize',24);title('PTT_P vs. Vessel Type');
subplot(224)
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
%   plot(tAll(2,round(OrgPTT_P(j)),j), PAll(2,round(OrgPTT_P(j)),j), 'ro');
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
%   plot(tAll(2,round(OrgPTT_U(j)),j), UAll(2,round(OrgPTT_U(j)),j), 'ro');
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
  %     plot(tAll(1,t_plot,j), PAll(1,t_plot,j), 'r', 'LineWidth', 2);
  plot(tAll(2,t_plot,j), PAll(2,t_plot,j), 'b', 'LineWidth', 2);
  %     plot(tAll(3,t_plot,j), PAll(3,t_plot,j), 'k', 'LineWidth', 2);
  ylabel('Area(um2)');xlabel('t(s)');
  %     legend('StartPt', 'MidPt', 'EndPt');
  title(['Area Vessel ' num2str(SegName(j),'%04d')]);
  
  saveas(h,['Vessel ' num2str(SegName(j),'%04d') '_' num2str(j,'%03d')]);
  saveas(h,['Vessel ' num2str(SegName(j),'%04d') '_' num2str(j,'%03d')],'png');
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
    subplot(3,2,i);hold on;
    plot(PAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
    plot(PpAll,'r','LineWidth',2);
    plot(PnAll,'k','LineWidth',2);
    subplot(3,2,i+3);hold on;
    plot(UAll(WIA_pt,t_plot,VesID),'b','LineWidth',2);
    plot(UpAll,'r','LineWidth',2);
    plot(UnAll,'k','LineWidth',2);
  end
%   saveas(h,['VesselWIA ' num2str(SegName(j),'%04d') '_' num2str(j,'%03d')],'png');
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
% MeanU(MeanU(:,2)<0,2)=-MeanU(MeanU(:,2)<0,2);
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

% ErrQ=sqrt(sum(((MeanQ(:,2)-FuncPara(:,8))./FuncPara(:,8)).^2)./VesNum);
% ErrP=sqrt(sum(((MeanP(:,2)-FuncPara(:,4))./FuncPara(:,4)).^2)./VesNum);
% ErrQ=mean(abs((MeanQ(:,2)-FuncPara(:,8))./FuncPara(:,8)));
% ErrP=mean(abs((MeanP(:,2)-FuncPara(:,4))./FuncPara(:,4)));
ErrQ=sum(abs((MeanQ(:,2)-FuncPara(:,8))./FuncPara(:,8)))/VesNum;
ErrP=sum(abs((MeanP(:,2)-FuncPara(:,4))./FuncPara(:,4)))/VesNum;
ErrA=mean(abs(sqrt(4/pi*(MeanA(:,2))-FuncPara(:,3))./FuncPara(:,3)));
ErrU=mean(abs((MeanU(:,2)-FuncPara(:,7))./FuncPara(:,7)));

%% Compare waveforms of different vessels
figure;hold on;
Ind=[1 2001];
tmpInd1=find(SegName==Ind(1));
tmpInd2=find(SegName==Ind(2));
if max(tmpInd1)>546
  fprintf(stderr,'Maximum vessel number is 546\n');
else
  % Scale first
  Tmp1=UAll(1,all_plot,tmpInd1);
  Tmp2=UAll(1,all_plot,tmpInd2);
  % Plot scaled curves
  plot(tAll(1,all_plot,tmpInd1),Tmp1,'.-');
  plot(tAll(1,all_plot,tmpInd2),Tmp2,'r.-');
end

if saveFlag
  clear saveFlag;
  save('546.mat');
end