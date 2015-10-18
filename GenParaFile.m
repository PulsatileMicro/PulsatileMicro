%% 生成符合Graph软件的parameter文件
% 格式为:
function GenParaFile(Para,SegName,ParaName,MinValue,MaxValue)
% 将Para转换为RGB颜色
colorOrder=16;
if nargin==3
  MapValue=uint8((Para-min(Para))./(max(Para)-min(Para))*colorOrder)+1;
elseif nargin==5
  MapValue=uint8((Para-MinValue)./(MaxValue-MinValue)*colorOrder)+1;
else
  error('Number of input arguments is wrong!\n');
end
MapValue(MapValue>colorOrder)=colorOrder;
cmap=colormap(jet(colorOrder));

% 将SegName, Para和RGB分量写入文件
fid=fopen([ParaName '.pra'],'w');
for i=1:length(Para)
  fprintf(fid,'%d\t%f\t%f\t%f\t%f\n',SegName(i),Para(i),cmap(MapValue(i),1),cmap(MapValue(i),2),cmap(MapValue(i),3));
end
fclose(fid);
end

