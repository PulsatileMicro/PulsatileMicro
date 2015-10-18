%% 生成1D模型输入文件(.in)的主函数
% clear;clc;close all;warning off;
% 
% ArtFile='I:\Kassab_Data\Set1\LAD\Set1_LAD_ART.txt';
% fid=fopen(ArtFile,'r');
% C=textscan(fid,'%d %d %f %f %f %f %f %f %f %f %d %d %d %d %d %d %d');
% fclose(fid);
% 1-输入类型, 0:q 1:p
ModelParam=[1];

VesNum=length(C{1});
SegName=C{1};
From=C{16};
To=C{17};
Len=C{3};
Diam=C{4};
Visc=2e-3*ones(VesNum,1);

%% 1D模型边界参数存储矩阵初始化
% 此处边界指的是所有血管段的边界，而非仅仅指血管网络的出入边界
% 所有血管段入口与出口处的边界类型，大小为VesNum*2
BCTypeAll=[];
% 对于分叉、汇聚、接合边界，Bifur为其子血管编号
% 对于输入边界，Bifur为其输入类型编号，如u 0, u 1的0和1
% 对于输出边界，Bifur不起作用
Bifur=zeros(VesNum,4);
% 对于分叉、汇聚、接合边界，BCVal不起作用
% 对于输入边界，BCVal为输入边界的值，如a=PI*1e-4
% 对于输出边界，BCVal为输出边界R或T的值
BCVal=zeros(VesNum,1);
% 边界流速矩阵
BoundData=[];

%% 为所有血管设置边界
%%%%%%%%%%%%%%%%%%%%%
%%%% 真实血管网络 %%%%
%%%%%%%%%%%%%%%%%%%%%
for i=1:VesNum
  BCType=[];  % 临时向量，存储边界类型
  %%% 分析输入边界
  inInd1=find(From(i)==From); % 查询其他从该段血管起点出发的血管
  inInd2=find(From(i)==To);   % 查询流入该段血管起点的血管
  if length(inInd1)==2
    % 如果有2条血管具有该起点，说明该节点是一个分叉节点
    % 此时必有且仅有一条血管流入该点，即inInd2
    BCType=[BCType 'B'];
    inInd1(inInd1==i)=[]; % 删除本身这条血管的编号
    Bifur(i,1:2)=[inInd2 inInd1];
  elseif length(inInd2)==2
    % 如果有2条血管流入该起点，说明该节点是一个汇聚节点
    BCType=[BCType 'C'];
    Bifur(i,1:2)=inInd2;
  elseif length(inInd2)==0
    % 如果没有血管流入该起点，说明该段血管是一个流入边界
    switch ModelParam(1)
      case 1
        BCType=[BCType 'q'];
        BCVal(i)=40000/60/1e12;
        Bifur(i,1:2)=[0.45 2];
      case 2
        BCType=[BCType 'p'];
        BCVal(i)=SS_Press(i);
        Bifur(i,1:2)=[0.45 2];
    end
    BoundData=[BoundData;BCVal(i) i 0];
  else
    % 以上情况均不符，则为接合血管
    BCType=[BCType 'J'];
    Bifur(i,1:2)=inInd2;
  end
  
  %%% 分析输出边界
  outInd1=find(To(i)==From);  % 查询该段血管流向的血管
  outInd2=find(To(i)==To);    % 查询其他到达该段血管终点的血管
  if length(outInd1)==2
    % 如果该段血管流向两根血管，说明该节点是一个分叉节点
    BCType=[BCType 'B'];
    Bifur(i,3:4)=outInd1;
  elseif length(outInd2)==2
    % 如果有2条血管到达该终点，说明该节点是一个汇聚节点
    BCType=[BCType 'C'];
    outInd2(outInd2==i)=[];
    Bifur(i,3:4)=[outInd2 outInd1];
  elseif length(outInd1)==0
    % 如果该终点不是任何血管的起点，说明该节点是一个流出边界
    % 稳态模型
    BCType=[BCType 'p'];
    BCVal(i)=SS_Press(i);
%   else
%     % 如果边界类型为流量，则设置流量为边界，符号为负
%     BCType=[BCType 'q'];
%     BCVal(i)=-40000/60/1e12;
  else
    % 以上情况均不符，则为接合血管
    BCType=[BCType 'J'];
    Bifur(i,3:4)=outInd1;
  end
  BCTypeAll=[BCTypeAll;BCType];
end