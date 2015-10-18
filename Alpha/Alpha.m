%%仿真数据的Order分析
clear all;
DatFile='T2810_h.dat';
PrnFile='T2810_h.prn';
[DataArray Boundary]=ReadData(DatFile);
VelNum=length(DataArray(:,1));
From=DataArray(:,3);
To=DataArray(:,4);
Diam=DataArray(:,5);

count=0;
for i=1:VelNum
  BifurIndex=find(From==To(i));
  if length(BifurIndex)==2
    count=count+1;
    if Diam(BifurIndex(1))<Diam(BifurIndex(2))
      Gamma=Diam(BifurIndex(1))/Diam(BifurIndex(2));
    else
      Gamma=Diam(BifurIndex(2))/Diam(BifurIndex(1));
    end
    Debug(count,:)=[i BifurIndex(1) BifurIndex(2) (Diam(BifurIndex(1))^2+Diam(BifurIndex(2))^2)/Diam(i)^2 Gamma];
  end
end
