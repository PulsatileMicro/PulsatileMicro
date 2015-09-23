% Determine the order of each vessel
% Order is defined as the number of bifurcations 
% to the main feeding arteriole
% clear;clc;close all;
DatFile='T2810_h.DAT';
PrnFile='T2810_h_adapted.prn';
[DataArray Boundary FuncPara]=ReadData(DatFile,PrnFile);
SegName=DataArray(:,1);
From=DataArray(:,3);
To=DataArray(:,4);
VesNum=length(SegName);
AllOrder=zeros(VesNum,1);
% This loop determines the order of art & cap
for i=1:VesNum
    tmpInd=i;
    order=1;
    while 1
        if From(tmpInd)==830;
            AllOrder(i)=order;
            break;
        else
            inInd=find(From(tmpInd)==To);
            if length(inInd)==1
                tmpInd=inInd;
                order=order+1;
            elseif length(inInd)==2
                order=0;
                break;
            else
                order=0;
                break;
            end
        end
    end
end
AllOrder(AllOrder==0)=max(AllOrder)+1;
AllOrder=max(AllOrder)-AllOrder+1;
% PlotMorph(AllOrder);
Order_RGB=ColorScal(AllOrder,256);
