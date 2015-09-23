%%作图test
clear all;clc;close all;
DatFile='T2810_h.dat';
PrnFile='T2810_h.prn';
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
VelType=DataArray(:,2);
VelNum=length(DataArray(:,1));
load('1D_Normal.mat');

% Diam=DataArray(:,5);
% DelInd=find(DataArray(:,1)>200 & DataArray(:,1)<300);
DelInd=[];
P=MeanP(:,2);
U=MeanU(:,2);
Uexp=DataArray(:,8);
Acnt=0;
Ccnt=0;
Vcnt=0;
P(DelInd)=[];U(DelInd)=[];VelType(DelInd)=[];Diam(DelInd)=[];
for i=1:VelNum-length(DelInd)
  if VelType(i)==1
    Acnt=Acnt+1;
    Art(Acnt,:)=[Diam(i) P(i) U(i) Uexp(i)];
  elseif VelType(i)==2
    Ccnt=Ccnt+1;
    Cap(Ccnt,:)=[Diam(i) P(i) U(i) Uexp(i)];
  else
    Vcnt=Vcnt+1;
    Vein(Vcnt,:)=[Diam(i) P(i) U(i) Uexp(i)];
  end
end
%Diam
Amax=max(Art(:,1));
Amin=min(Art(:,1));
Cmax=max(Cap(:,1));
Cmin=min(Cap(:,1));
Vemax=max(Vein(:,1));
Vemin=min(Vein(:,1));

%%Art
Spa=6;  %管径间隔
AOrder=floor(Amax/Spa);
for i=1:AOrder
  Index1=find(Art(:,1)<8+(AOrder-i)*Spa+3);
  Index2=find(Art(:,1)>=8+(AOrder-i-1)*Spa+3);
  Index=intersect(Index1,Index2);
  if isempty(Index)
    continue;
  end
  %     8+(AOrder-i)*Spa+3
  %     8+(AOrder-i-1)*Spa+3
  Pmean=mean(Art(Index,2));
  Pmax=max(Art(Index,2));
  Pmin=min(Art(Index,2));
  Pstd=std(Art(Index,2));
  Vmean=mean(Art(Index,3));
  Vmax=max(Art(Index,3));
  Vmin=min(Art(Index,3));
  Vstd=std(Art(Index,3));
  Vexp_mean=mean(Art(Index,4));
  Vexp_max=max(Art(Index,4));
  Vexp_min=min(Art(Index,4));
  Vexp_std=std(Art(Index,4));
  APPoint(i,:)=[Pmin Pmean Pmax Pstd];
  AVelPoint(i,:)=[Vmin Vmean Vmax Vstd];
  AVelexpPoint(i,:)=[Vexp_min Vexp_mean Vexp_max Vexp_std];
end

%%Vein
VOrder=floor(Vemax/Spa)-2;
for i=1:VOrder
  Index1=find(Vein(:,1)<8+(i-1)*Spa+3);
  Index2=find(Vein(:,1)>=8+(i-2)*Spa+3);
  Index=intersect(Index1,Index2);
  8+(i-1)*Spa+3
  8+(i-2)*Spa+3
  Pmean=mean(Vein(Index,2));
  Pmax=max(Vein(Index,2));
  Pmin=min(Vein(Index,2));
  Pstd=std(Vein(Index,2));
  Vmean=mean(Vein(Index,3));
  Vmax=max(Vein(Index,3));
  Vmin=min(Vein(Index,3));
  Vstd=std(Vein(Index,3));
  Vexp_mean=mean(Vein(Index,4));
  Vexp_max=max(Vein(Index,4));
  Vexp_min=min(Vein(Index,4));
  Vexp_std=std(Vein(Index,4));
  VPPoint(i,:)=[Pmin Pmean Pmax Pstd];
  VVelPoint(i,:)=[Vmin Vmean Vmax Vstd];
  VVelexpPoint(i,:)=[Vexp_min Vexp_mean Vexp_max Vexp_std];
end

% 双坐标轴的版本
% x=1:AOrder+VOrder;
% PPoint=[APPoint;VPPoint];
% VelPoint=[AVelPoint;VVelPoint];
% VelexpPoint=[AVelexpPoint;VVelexpPoint];
% [AX,H1,H2]=plotyy(x,PPoint(:,2),x,VelPoint(:,2),'plot');
% set(get(AX(1),'Ylabel'),'String','Pressure(mmHg)')
% set(get(AX(2),'Ylabel'),'String','Velocity(mm/s)')
% xlabel('Diameter(um)');
% set(H1,'linestyle','-','color','k');
% set(H2,'linestyle','- -','color','k');
% legend([H1 H2],'Pressure','Velocity') %标注两条线图例
% legend('boxoff')%去掉图例边框
% box off;%去掉图边框
% set(AX(:),'Ycolor','k') %设定两个Y轴的颜色为黑色
% set(AX(1),'YLim',[0 100],'ytick',[0:20:100]); %设置y轴间隔
% % set(AX(1),'YLim',[0 12],'ytick',[0:3:12]);
% set(AX(2),'YLim',[0 12],'ytick',[0:3:12]);
% set(AX,'XTick',[-3:1:AOrder+VOrder],'XLim',[-3 AOrder+VOrder],'XLimMode','manual','xTicklabel',{'50','44','38','32','26','20','14','8','8','14','20','26','32','38','44','50'}) % 设置x轴范围

% 单坐标轴、子图的版本
x=1:AOrder+VOrder;
PPoint=[APPoint;VPPoint];
VelPoint=[AVelPoint;VVelPoint];
VelexpPoint=[AVelexpPoint;VVelexpPoint];
figure;
subplot(121);
errorbar(x,PPoint(:,2),PPoint(:,4),'k.-','LineWidth',2);
xlabel('Vessel Diameter(\mum)');
ylabel('Pressure(mmHg)');
subplot(122);
errorbar(x,VelPoint(:,2),VelPoint(:,4),'k.-','LineWidth',2);
xlabel('Vessel Diameter(\mum)');
ylabel('Velocity(mm/s)');