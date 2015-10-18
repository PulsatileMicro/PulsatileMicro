%% ���ɷ���Graph�����parameter�ļ�
% ��ʽΪ:
function GenParaFile(Para,SegName,ParaName,MinValue,MaxValue)
% ��Paraת��ΪRGB��ɫ
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

% ��SegName, Para��RGB����д���ļ�
fid=fopen([ParaName '.pra'],'w');
for i=1:length(Para)
  fprintf(fid,'%d\t%f\t%f\t%f\t%f\n',SegName(i),Para(i),cmap(MapValue(i),1),cmap(MapValue(i),2),cmap(MapValue(i),3));
end
fclose(fid);
end

