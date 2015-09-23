%%仿真数据的Order分析
clear all;
DatFile='T2810_h.dat';
PrnFile='T2810_h.prn';
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
% [FuncPara]=textread('T_2810.txt');
VelNum=length(DataArray(:,1));
From=DataArray(:,3);
To=DataArray(:,4);
DiamInitial=FuncPara(:,2);
DiamAdapted=FuncPara(:,3);
Flow=FuncPara(:,8);
From1=From;
To1=To;
Index=find(Flow<0);
From(Index)=To1(Index);
To(Index)=From1(Index);

count=0;
for i=1:VelNum
    BifurIndex=find(From==To(i));
    if length(BifurIndex)==2
        count=count+1;
%         Debug(count,:)=[i BifurIndex(1) BifurIndex(2) (Diam(BifurIndex(1))^2+Diam(BifurIndex(2))^2)/Diam(i)^2];
        if DiamInitial(BifurIndex(1))>DiamInitial(BifurIndex(2))
            DebugInitial(count,:)=[DiamInitial(BifurIndex(2))/DiamInitial(BifurIndex(1)) DiamInitial(BifurIndex(2))/DiamInitial(i) DiamInitial(BifurIndex(1))/DiamInitial(i)];
        else
            DebugInitial(count,:)=[DiamInitial(BifurIndex(1))/DiamInitial(BifurIndex(2)) DiamInitial(BifurIndex(1))/DiamInitial(i) DiamInitial(BifurIndex(2))/DiamInitial(i)];
        end
        if DiamAdapted(BifurIndex(1))>DiamAdapted(BifurIndex(2))
            DebugAdapted(count,:)=[DiamAdapted(BifurIndex(2))/DiamAdapted(BifurIndex(1)) DiamAdapted(BifurIndex(2))/DiamAdapted(i) DiamAdapted(BifurIndex(1))/DiamAdapted(i)];
        else
            DebugAdapted(count,:)=[DiamAdapted(BifurIndex(1))/DiamAdapted(BifurIndex(2)) DiamAdapted(BifurIndex(1))/DiamAdapted(i) DiamAdapted(BifurIndex(2))/DiamAdapted(i)];
        end
    end
end

x=0:0.01:1;
y1=(1+x.^3).^(-1/3);
y2=x.*(1+x.^3).^(-1/3);

[a,b]=sort(DebugInitial(:,1),'ascend');
Debug1=DebugInitial(b,:);
pI3=polyfit(Debug1(:,1),Debug1(:,3),3);
pI2=polyfit(Debug1(:,1),Debug1(:,2),3);
xcurve=0:0.05:1;
ycurveI3=polyval(pI3,xcurve);
ycurveI2=polyval(pI2,xcurve);
figure;
subplot(2,1,1);
plot(DebugInitial(:,1),DebugInitial(:,3),'ko');hold on
plot(x,y1,'k-');hold on
plot(xcurve,ycurveI3,'r-');
xlabel('D2/D1');
ylabel('D1/D0');
title('Initial546');

subplot(2,1,2);
plot(DebugInitial(:,1),DebugInitial(:,2),'ko');hold on
plot(x,y2,'k-');hold on
plot(xcurve,ycurveI2,'r-');
xlabel('D2/D1');
ylabel('D2/D0');


[a,b]=sort(DebugAdapted(:,1),'ascend');
Debug2=DebugAdapted(b,:);
pA3=polyfit(Debug2(:,1),Debug2(:,3),3);
pA2=polyfit(Debug2(:,1),Debug2(:,2),3);
xcurve=0:0.05:1;
ycurveA3=polyval(pA3,xcurve);
ycurveA2=polyval(pA2,xcurve);
figure;
subplot(2,1,1);
plot(DebugAdapted(:,1),DebugAdapted(:,3),'ko');hold on
plot(x,y1,'k-');hold on
plot(xcurve,ycurveA3,'r-');
xlabel('D2/D1');
ylabel('D1/D0');
title('Adapted546');

subplot(2,1,2);
plot(DebugAdapted(:,1),DebugAdapted(:,2),'ko');hold on
plot(x,y2,'k-');hold on
plot(xcurve,ycurveA2,'r-');
xlabel('D2/D1');
ylabel('D2/D0');


DatFile='T15_8V_h.dat';
[DataArray2 Boundary2]=ReadData1(DatFile);
VelNum2=length(DataArray2(:,1));
From2=DataArray2(:,3);
To2=DataArray2(:,4);
DiamInitial2=DataArray2(:,5);
count2=0;
for i=1:VelNum2
    BifurIndex=find(From2==To2(i));
    if length(BifurIndex)==2
        count2=count2+1;
        if DiamInitial2(BifurIndex(1))>DiamInitial2(BifurIndex(2))
            DebugInitial2(count2,:)=[DiamInitial2(BifurIndex(2))/DiamInitial2(BifurIndex(1)) DiamInitial2(BifurIndex(2))/DiamInitial2(i) DiamInitial2(BifurIndex(1))/DiamInitial2(i)];
        else
            DebugInitial2(count2,:)=[DiamInitial2(BifurIndex(1))/DiamInitial2(BifurIndex(2)) DiamInitial2(BifurIndex(1))/DiamInitial2(i) DiamInitial2(BifurIndex(2))/DiamInitial2(i)];
        end
    end
end

DatFile='T1_8o_h.dat';    %边界节点非连续编号，不能使用边界节点重排模块
[DataArray3 Boundary3]=ReadData1(DatFile);
VelNum3=length(DataArray3(:,1));
From3=DataArray3(:,3);
To3=DataArray3(:,4);
DiamInitial3=DataArray3(:,5);
count3=0;
for i=1:VelNum3
    BifurIndex=find(From3==To3(i));
    if length(BifurIndex)==2
        count3=count3+1;
        if DiamInitial3(BifurIndex(1))>DiamInitial3(BifurIndex(2))
            DebugInitial3(count3,:)=[DiamInitial3(BifurIndex(2))/DiamInitial3(BifurIndex(1)) DiamInitial3(BifurIndex(2))/DiamInitial3(i) DiamInitial3(BifurIndex(1))/DiamInitial3(i)];
        else
            DebugInitial3(count3,:)=[DiamInitial3(BifurIndex(1))/DiamInitial3(BifurIndex(2)) DiamInitial3(BifurIndex(1))/DiamInitial3(i) DiamInitial3(BifurIndex(2))/DiamInitial3(i)];
        end
    end
end

DebugAll=[DebugInitial; DebugInitial2; DebugInitial3];
x=0:0.01:1;
y1=(1+x.^3).^(-1/3);
y2=x.*(1+x.^3).^(-1/3);

[a,b]=sort(DebugAll(:,1),'ascend');
Debug3=DebugAll(b,:);
pAl3=polyfit(Debug3(:,1),Debug3(:,3),3);
pAl2=polyfit(Debug3(:,1),Debug3(:,2),3);
xcurve=0:0.05:1;
ycurveAl3=polyval(pAl3,xcurve);
ycurveAl2=polyval(pAl2,xcurve);
figure;
subplot(2,1,1);
plot(DebugAll(:,1),DebugAll(:,3),'ko');hold on
plot(x,y1,'k-');hold on
plot(xcurve,ycurveAl3,'r-');
xlabel('D2/D1');
ylabel('D1/D0');
title('All Networks');

subplot(2,1,2);
plot(DebugAll(:,1),DebugAll(:,2),'ko');hold on
plot(x,y2,'k-');hold on
plot(xcurve,ycurveAl2,'r-');
xlabel('D2/D1');
ylabel('D2/D0');
