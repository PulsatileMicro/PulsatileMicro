%% Steady State模型仿真结果分析
% Before running this script, what you have to do are:
% 1. Run main_1D.m (in the directory with main_1D.m file)
% 2. Enter the data directory (such as 546_qR_E1)
% 3. Run this script and click "add to path" but not "change directory"
close all;clc;
QAll=zeros(VesNum,1);
PAll=zeros(VesNum,1);
for j=1:VesNum
  if VesNum == 1
    fileName = [NetTypeName '.his'];
  else
    fileName = [NetTypeName '_' int2str(j) '.his'];
  end
  
  fid=fopen(fileName, 'r');
  C=textscan(fid, '%f %f', 'HeaderLines',2);
  QAll(j)=C{1};
  PAll(j)=C{2};
  fclose(fid);
end