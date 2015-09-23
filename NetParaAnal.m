%% 血管网络结构参数分析
clear;clc;close all;
% 参与分析的网络
% 1-Net_546网络
% 2-Egg_636网络
% 3-CAM网络

%% Net_546
DatFile='T2810_h.DAT';
PrnFile='';
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
% 管径
Diam=DataArray(:,5);
% 长度
Len=DataArray(:,6);
% 壁厚
VesType=DataArray(:,2);
VesNum=length(Diam);
WallTh=zeros(VesNum,1);
From=DataArray(:,3);
To=DataArray(:,4);
for i=1:VesNum
  if VesType(i)==1
    WallTh(i)=Diam(i)*0.2662;
  elseif VesType(i)==2
    WallTh(i)=Diam(i)*0.1812;
  elseif VesType(i)==3
    WallTh(i)=Diam(i)*0.0686;
  else
    WallTh(i)=Diam(i)*0.2;
  end
end
% 杨氏模量
E=zeros(VesNum,1);
for i=1:VesNum
  if VesType(i)==1
    E(i)=3.5e5;
  elseif VesType(i)==3
    E(i)=3.88e5;
  else
    E(i)=3.7e5;
  end
end
% 粘滞度
load 546_Meas.mat;

% 计算各类参数
VenInd_546=(VesType==3);
ArtCapInd_546=(VesType~=3);
Res_546=CalcVesRes(Diam,Visc,Len);
f=E.*WallTh.*Diam/2./(24.*pi.*Visc.*Len);
LogInvF_546=log10(1./f);
% LogInvF_546(VesType==3)=min(LogInvF_546);
EHPIETA_546=log10(E.*WallTh.*Visc*24);
% EHPIETA_546(VenInd_546)=min(EHPIETA_546);
LRRatio_546=log10(Len./Diam);
% LRRatio_546(VenInd_546)=min(LRRatio_546);
CRRatio_546=log10(3*pi^2*(Diam/2).^7./(16*Visc.*Len.*E.*WallTh));
% CRRatio_546(VenInd_546)=min(CRRatio_546);
[VesOrder OrderRange]=CalcVesOrder(From,To,830);
for i=1:length(OrderRange)-1
  ind=(VesOrder==i);
  Debug(i,1)=mean(LRRatio_546(ind));
end
figure;plot(Debug);

%% Egg_636
DatFile='Egg636_morph_update.dat';
PrnFile='';
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
% 管径
Diam=DataArray(:,5);
% 长度
Len=DataArray(:,6);
% 壁厚
VesType=DataArray(:,2);
VesNum=length(Diam);
WallTh=zeros(VesNum,1);
From=DataArray(:,3);
To=DataArray(:,4);
for i=1:VesNum
  if VesType(i)==1
    WallTh(i)=Diam(i)*0.2662;
  elseif VesType(i)==2
    WallTh(i)=Diam(i)*0.1812;
  elseif VesType(i)==3
    WallTh(i)=Diam(i)*0.0686;
  else
    WallTh(i)=Diam(i)*0.2;
  end
end
% 杨氏模量
E=zeros(VesNum,1);
for i=1:VesNum
  if VesType(i)==1
    E(i)=3.5e5;
  elseif VesType(i)==3
    E(i)=3.88e5;
  else
    E(i)=3.7e5;
  end
end
% 粘滞度
load Egg_636_Meas.mat;

% 计算各类参数
VenInd_636=(VesType==3);
ArtCapInd_636=(VesType~=3);
Res_636=CalcVesRes(Diam,Visc,Len);
f=E.*WallTh.*Diam/2./(24.*pi.*Visc.*Len);
LogInvF_636=log10(1./f);
% LogInvF_636(VesType==3)=min(LogInvF_636);
EHPIETA_636=log10(E.*WallTh.*Visc*24);
% EHPIETA_636(VenInd_636)=min(EHPIETA_636);
LRRatio_636=log10(Len./Diam);
% LRRatio_636(VenInd_636)=min(LRRatio_636);
CRRatio_636=log10(3*pi^2*(Diam/2).^7./(16*Visc.*Len.*E.*WallTh));
% CRRatio_636(VenInd_636)=min(CRRatio_636);
[VesOrder OrderRange]=CalcVesOrder(From,To,137);
for i=1:length(OrderRange)-1
  ind=(VesOrder==i);
  Debug(i,1)=mean(f(ind));
end
figure;plot(Debug);

%% Sub_CAM
DatFile='SubCAM.DAT';
PrnFile='';
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
% 管径
Diam=DataArray(:,5);
% 长度
Len=DataArray(:,6);
% 壁厚
VesType=DataArray(:,2);
VesNum=length(Diam);
WallTh=zeros(VesNum,1);
From=DataArray(:,3);
To=DataArray(:,4);
for i=1:VesNum
  if VesType(i)==1
    WallTh(i)=Diam(i)*0.2662;
  elseif VesType(i)==2
    WallTh(i)=Diam(i)*0.1812;
  elseif VesType(i)==3
    WallTh(i)=Diam(i)*0.0686;
  else
    WallTh(i)=Diam(i)*0.2;
  end
end
% 杨氏模量
E=zeros(VesNum,1);
for i=1:VesNum
  if VesType(i)==1
    E(i)=3.5e5;
  elseif VesType(i)==3
    E(i)=3.88e5;
  else
    E(i)=3.7e5;
  end
end
% 粘滞度
load subCAM_Meas.mat;
% 计算各类参数
VenInd_CAM=(VesType==3);
ArtCapInd_CAM=(VesType~=3);
f=E.*WallTh.*Diam/2./(24.*pi.*Visc.*Len);
LogInvF_CAM=log10(1./f);
% LogInvF_CAM(VesType==3)=min(LogInvF_CAM);
EHPIETA_CAM=log10(E.*WallTh.*Visc*24);
% EHPIETA_CAM(VenInd_CAM)=min(EHPIETA_CAM);
LRRatio_CAM=log10(Len./Diam);
% LRRatio_CAM(VenInd_CAM)=min(LRRatio_CAM);
CRRatio_CAM=log10(3*pi^2*(Diam/2).^7./(16*Visc.*Len.*E.*WallTh));
% CRRatio_CAM(VenInd_CAM)=min(CRRatio_CAM);
[VesOrder OrderRange]=CalcVesOrder(From,To,4956);
for i=1:length(OrderRange)-1
  ind=(VesOrder==i);
  Debug(i,1)=mean(f(ind));
end
figure;plot(Debug);

MaxLogInvF=max([LogInvF_546;LogInvF_636;LogInvF_CAM]);
MinLogInvF=min([LogInvF_546;LogInvF_636;LogInvF_CAM]);
MaxEHPIETA=max([EHPIETA_546;EHPIETA_636;EHPIETA_CAM]);
MinEHPIETA=min([EHPIETA_546;EHPIETA_636;EHPIETA_CAM]);
MaxLRRatio=max([LRRatio_546;LRRatio_636;LRRatio_CAM]);
MinLRRatio=min([LRRatio_546;LRRatio_636;LRRatio_CAM]);
MaxCRRatio=max([CRRatio_546;CRRatio_636;CRRatio_CAM]);
MinCRRatio=min([CRRatio_546;CRRatio_636;CRRatio_CAM]);

% LogInvF_546(VenInd_546)=MinLogInvF;
% LogInvF_636(VenInd_636)=MinLogInvF;
% LogInvF_CAM(VenInd_CAM)=MinLogInvF;
% EHPIETA_546(VenInd_546)=MinEHPIETA;
% EHPIETA_636(VenInd_636)=MinEHPIETA;
% EHPIETA_CAM(VenInd_CAM)=MinEHPIETA;
% LRRatio_546(VenInd_546)=MinLRRatio;
% LRRatio_636(VenInd_636)=MinLRRatio;
% LRRatio_CAM(VenInd_CAM)=MinLRRatio;
% CRRatio_546(VenInd_546)=MinCRRatio;
% CRRatio_636(VenInd_636)=MinCRRatio;
% CRRatio_CAM(VenInd_CAM)=MinCRRatio;

% PlotTopol(LogInvF_546,1,MinLogInvF,MaxLogInvF);
% PlotTopol(LogInvF_636,9,MinLogInvF,MaxLogInvF);
% PlotTopol(LogInvF_CAM,8,MinLogInvF,MaxLogInvF);
% PlotTopol(EHPIETA_546,1,MinEHPIETA,MaxEHPIETA);
% PlotTopol(EHPIETA_636,9,MinEHPIETA,MaxEHPIETA);
% PlotTopol(EHPIETA_CAM,8,MinEHPIETA,MaxEHPIETA);
PlotTopol(LRRatio_546,1,MinLRRatio,MaxLRRatio);
PlotTopol(LRRatio_636,9,MinLRRatio,MaxLRRatio);
PlotTopol(LRRatio_CAM,8,MinLRRatio,MaxLRRatio);
% PlotTopol(CRRatio_546,1,MinCRRatio,MaxCRRatio);
% PlotTopol(CRRatio_636,9,MinCRRatio,MaxCRRatio);
% PlotTopol(CRRatio_CAM,8,MinCRRatio,MaxCRRatio);
