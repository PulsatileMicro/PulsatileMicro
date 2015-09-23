% 画图程序：网络拓扑图 （参数分布） %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 血管网络类型 %%%%
function PlotTopol(Input,VesType,Title,MinValue,MaxValue)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID

switch VesType
  case {Net_546_ID, Net_546_Meas_ID}
    %546 network
    topolData='Men_546.opt';
    ElemNum=1847;
  case Egg_818_ID
    %Egg_818 network
    topolData='Vit_818.opt';
    ElemNum=4213;
  case Net_122_ID
    %122 network
    error(' 122 network has no topol data file !');
  case Net_389_ID
    %389 network
    topolData='Men_389.opt';
    ElemNum=3573;
  case Net_913_ID
    %913 network
    topolData='Men_913.opt';
    ElemNum=4509;
  case Egg_CAM_ID
    %CAM network
    topolData='CAM_7128.opt';
    ElemNum=18480;
  case Sub_CAM_ID
    topolData='SubCAM_opt.txt';
    ElemNum=2849;
  case Egg_636_ID
    topolData='Vit_636.opt';
    ElemNum=2529;
  otherwise
    error('Wrong setting of "VesType" !');
end

%读入topology数据
fid=fopen(topolData);
if VesType~=Net_389_ID
  tData=textscan(fid,'%d %d %s %d %d %d %d %d %d %d %d %d %f %f %f %f %d',ElemNum,'delimiter',',','HeaderLines',1);
else
  tData=textscan(fid,'%d %d %s %d %d %d %d %d %d %d %d %d',ElemNum,'delimiter','\t','HeaderLines',1);
end
SegmentId=tData{1};
ProxNodeX=tData{6};
ProxNodeY=tData{7};
DistNodeX=tData{10};
DistNodeY=tData{11};
if VesType~=Net_389_ID
  Diameter=tData{13};
else
  Diameter=30*ones(ElemNum,1);
end

colorOrder=256;
if nargin==3
  MapValue=uint8((Input-min(Input))./(max(Input)-min(Input))*colorOrder)+1;
elseif nargin==5
  MapValue=uint8((Input-MinValue)./(MaxValue-MinValue)*colorOrder)+1;
else
  error('Number of input arguments is wrong!\n');
end
MapValue(MapValue>colorOrder)=colorOrder;

figure;hold on;
whitebg('k');
cmap=colormap(jet(colorOrder));
%计数（通过计数建立Input与Elem的关系）
%Input对应完整血管段，Elem对应离散血管段
SegCnt=1;
if VesType==Egg_CAM_ID
  LastSegId=4;
elseif VesType==Egg_818_ID
  LastSegId=3;
elseif VesType==Sub_CAM_ID
%   LastSegId=1124;
  LastSegId=7422;
else
  LastSegId=1;
end

for i=1:ElemNum
  if LastSegId~=SegmentId(i)
    SegCnt=SegCnt+1;
    LastSegId=SegmentId(i);
  end
  tmpX=[ProxNodeX(i) DistNodeX(i)];
  tmpY=[-ProxNodeY(i) -DistNodeY(i)];
  plot(tmpX,tmpY,'LineWidth',Diameter(i)/15,'Color',cmap(MapValue(SegCnt),:),'Marker','o','MarkerSize',Diameter(i)/45);
end
title(Title);
