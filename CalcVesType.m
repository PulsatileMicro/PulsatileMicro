%%Classification for network vessels:arteries venules capillaries
% clear all;clc;
% Network data file name
DatFile='T2810_h.DAT';
PrnFile='T2810_h_adapted.prn';
% Read data from input files
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
VesNum=length(DataArray(:,1));
SegName=DataArray(:,1);
From=DataArray(:,3);
To=DataArray(:,4);
%修改某些血管方向
% Invert the output boundaries
InvInd=[2205 1215 215 216 1716 2275 1230 1231 231 1222 222 1235 1237 1238 238 237 2088 2082 2035];
for i=1:VesNum
  if (From(i)>800 && From(i)<900 && From(i)~=830) || ...
      SegName(i)==2205 || SegName(i)==1215 || ...
      SegName(i)==215 || SegName(i)==216 || SegName(i)==1716 ||...
      SegName(i)==2275 || SegName(i)==1230 || SegName(i)==1231 ||...
      SegName(i)==231 || SegName(i)==1222 || SegName(i)==222 ||...
      SegName(i)==1235 || SegName(i)==1237 || SegName(i)==1238 ||...
      SegName(i)==238 || SegName(i)==237 || SegName(i)==2088 || ...
      SegName(i)==2082 || SegName(i)==2035
    temp=From(i);
    From(i)=To(i);
    To(i)=temp;
    %        Flow(i)=-Flow(i);
  end
end
VesselType=zeros(VesNum,6);
for j=1:5
  for i=1:VesNum
    %动脉不直接与静脉相连
    %没有上一级入血管为动脉
    %没有下一级出血管为静脉
    %有两个上一级如血管为静脉
    %一个入血管，两个出血管，与入血管相同
    %一个入血管，一个出血管，入动脉，则为毛细血管；入静脉则为静脉
    BFIndex=find(Boundary(:,1)==From(i));
    BTIndex=find(Boundary(:,1)==To(i));  %Boundary vessels set to be capillaries;
    BFNum=length(BFIndex);
    BTNum=length(BTIndex);
    ToIndex=find(From==To(i));
    TToIndex=find(To==To(i));
    FromIndex=find(To==From(i));
    FFromIndex=find(From==From(i));
    ToNum=length(ToIndex);
    TToNum=length(TToIndex);
    FromNum=length(FromIndex);
    FFromNum=length(FFromIndex);
    VesselType(i,3)=ToNum;
    VesselType(i,4)=TToNum;
    VesselType(i,5)=FromNum;
    VesselType(i,6)=FFromNum;
    if BFNum==1
      VesselType(i,1)=1; %arteries
    elseif BTNum==1
      VesselType(i,1)=3; %
    elseif FromNum==2
      VesselType(i,1)=3;
    elseif FromNum==1&&ToNum==2
      if VesselType(FromIndex(1),1)==1&&VesselType(ToIndex(1),1)==3
        VesselType(i,1)=2;  %
      elseif  VesselType(FromIndex(1),1)==1&&VesselType(ToIndex(2),1)==3;  %
        VesselType(i,1)=2;  %
      elseif VesselType(FromIndex(1),1)==1
        VesselType(i,1)=1;  %
      elseif VesselType(FromIndex(1),1)==3
        VesselType(i,1)=3;
      else
        VesselType(i,1)=0;
      end
    elseif ToNum==1&&FromNum==1
      if VesselType(FromIndex(1),1)==1
        VesselType(i,1)=2; %capillaries
      elseif VesselType(FromIndex(1),1)==3
        VesselType(i,1)=3; %
      else
        VesselType(i,1)=0; %
      end
    else
      VesselType(i,1)=2; %arteries
    end
  end
end  %血管顺序，保证不出现Type=0情况
VesselType(:,2)=DataArray(:,2);

%人工调节不合适血管Type
Index=find(DataArray(:,1)==2047);
VesselType(Index,1)=3;
Index=find(DataArray(:,1)==2167);
VesselType(Index,1)=3;

%%画图
ElemNum=1847;
[SegmentId,...
  SegmentPartID,...
  SegmentName,...
  SegmentType,...
  ProxNodeNumber,...
  ProxNodeX,...
  ProxNodeY,...
  ProxNodeZ,...
  DistNodeNumber,...
  DistNodeX,...
  DistNodeY,...
  DistNodeZ,...
  Diameter,...
  Length,...
  Hematocrit,...
  Velocity,...
  ExpID]...
  =textread('T28102000_opt.txt','%d %d %s %d %d %d %d %d %d %d %d %d %f %f %f %f %d',ElemNum,'delimiter',',');

colorOrder=64;
% MapValue=uint8((Input-min(Input))./(max(Input)-min(Input))*colorOrder)+1;
% MapValue(MapValue>colorOrder)=colorOrder;

figure;hold on;
whitebg('k');
cmap=colormap(jet(colorOrder));
SegCnt=1;

for i=1:ElemNum
  tmpX=[ProxNodeX(i) DistNodeX(i)];
  tmpY=[-ProxNodeY(i) -DistNodeY(i)];
  %    plot(tmpX,tmpY,'LineWidth',Diameter(i)/8,'Color',cmap(MapValue(SegCnt),:),'Marker','o','MarkerSize',Diameter(i)/45);
  %     plot(tmpX,tmpY,'w','LineWidth',Diameter(i)/12);
  %     if FuncPara(SegCnt,19)==1
  %         plot(tmpX,tmpY,'r','LineWidth',Diameter(i)/8,'Marker','o','MarkerSize',Diameter(i)/40);
  %     elseif FuncPara(SegCnt,19)==3
  %         plot(tmpX,tmpY,'b','LineWidth',Diameter(i)/8,'Marker','o','MarkerSize',Diameter(i)/40);
  %     else
  %         plot(tmpX,tmpY,'y','LineWidth',Diameter(i)/8,'Marker','o','MarkerSize',Diameter(i)/40);
  %     end
  
  %     plot(tmpX,tmpY,'w','LineWidth',Diameter(i)/12);
  if VesselType(SegCnt,1)==1
    plot(tmpX,tmpY,'r','LineWidth',Diameter(i)/8,'Marker','o','MarkerSize',Diameter(i)/40);
  elseif VesselType(SegCnt,1)==3
    plot(tmpX,tmpY,'b','LineWidth',Diameter(i)/8,'Marker','o','MarkerSize',Diameter(i)/40);
  else
    plot(tmpX,tmpY,'y','LineWidth',Diameter(i)/8,'Marker','o','MarkerSize',Diameter(i)/40);
  end
  if SegmentPartID(i) == 99
    SegCnt=SegCnt+1;
  end
end