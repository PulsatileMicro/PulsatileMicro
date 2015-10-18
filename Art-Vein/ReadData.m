%% Read the .dat and .prn file
function [DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile)
DataArray=[];
Boundary=[];
FuncPara=[];
% Read topology information
fiddat=fopen(DatFile,'r');
% 	Topology=textscan(fiddat,'%u %u %u %u %f %f %f %f','HeaderLines',6);
for i=1:7
	line=fgetl(fiddat);
end
while line~=-1
	if length(line)<10
		break;
	else
		temp=strread(line);
		DataArray=[DataArray;temp];
		line=fgetl(fiddat);
	end
end
fgetl(fiddat);	% omit the headerlines
line=fgetl(fiddat);
while line~=-1
	temp=strread(line);
	Boundary=[Boundary;temp];
	line=fgetl(fiddat);
end
fclose(fiddat);

% Read functional information
fidprn=fopen(PrnFile,'r');
line=fgetl(fidprn);
for i=1:4
	line=fgetl(fidprn);
end
while line~=-1
	if length(line)<10
		break;
	else
		temp=strread(line);
		FuncPara=[FuncPara;temp];
		line=fgetl(fidprn);
	end
end
fclose(fidprn);
