%%仿真数据的Order分析
clear all;
DatFile='T2810_h.dat';
PrnFile='T2810_h_adapted.prn';
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
VelNum=length(DataArray(:,1));
From=DataArray(:,3);
To=DataArray(:,4);
Diam=FuncPara(:,3);

count=0;
for i=1:VelNum
  BifurIndex=find(From==To(i));
  if length(BifurIndex)==2
    count=count+1;
    Debug(count,:)=[i
      BifurIndex(1)
      BifurIndex(2)
      (Diam(BifurIndex(1))^2+Diam(BifurIndex(2))^2)/Diam(i)^2
      Diam(i)^3/(Diam(BifurIndex(1))^3+Diam(BifurIndex(2))^3)];
  end
end
