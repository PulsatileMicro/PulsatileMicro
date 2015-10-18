%% ����1Dģ�������ļ�(.in)��������
% clear;clc;close all;warning off;
% 
% ArtFile='I:\Kassab_Data\Set1\LAD\Set1_LAD_ART.txt';
% fid=fopen(ArtFile,'r');
% C=textscan(fid,'%d %d %f %f %f %f %f %f %f %f %d %d %d %d %d %d %d');
% fclose(fid);
% 1-��������, 0:q 1:p
ModelParam=[1];

VesNum=length(C{1});
SegName=C{1};
From=C{16};
To=C{17};
Len=C{3};
Diam=C{4};
Visc=2e-3*ones(VesNum,1);

%% 1Dģ�ͱ߽�����洢�����ʼ��
% �˴��߽�ָ��������Ѫ�ܶεı߽磬���ǽ���ָѪ������ĳ���߽�
% ����Ѫ�ܶ��������ڴ��ı߽����ͣ���СΪVesNum*2
BCTypeAll=[];
% ���ڷֲ桢��ۡ��Ӻϱ߽磬BifurΪ����Ѫ�ܱ��
% ��������߽磬BifurΪ���������ͱ�ţ���u 0, u 1��0��1
% ��������߽磬Bifur��������
Bifur=zeros(VesNum,4);
% ���ڷֲ桢��ۡ��Ӻϱ߽磬BCVal��������
% ��������߽磬BCValΪ����߽��ֵ����a=PI*1e-4
% ��������߽磬BCValΪ����߽�R��T��ֵ
BCVal=zeros(VesNum,1);
% �߽����پ���
BoundData=[];

%% Ϊ����Ѫ�����ñ߽�
%%%%%%%%%%%%%%%%%%%%%
%%%% ��ʵѪ������ %%%%
%%%%%%%%%%%%%%%%%%%%%
for i=1:VesNum
  BCType=[];  % ��ʱ�������洢�߽�����
  %%% ��������߽�
  inInd1=find(From(i)==From); % ��ѯ�����Ӹö�Ѫ����������Ѫ��
  inInd2=find(From(i)==To);   % ��ѯ����ö�Ѫ������Ѫ��
  if length(inInd1)==2
    % �����2��Ѫ�ܾ��и���㣬˵���ýڵ���һ���ֲ�ڵ�
    % ��ʱ�����ҽ���һ��Ѫ������õ㣬��inInd2
    BCType=[BCType 'B'];
    inInd1(inInd1==i)=[]; % ɾ����������Ѫ�ܵı��
    Bifur(i,1:2)=[inInd2 inInd1];
  elseif length(inInd2)==2
    % �����2��Ѫ���������㣬˵���ýڵ���һ����۽ڵ�
    BCType=[BCType 'C'];
    Bifur(i,1:2)=inInd2;
  elseif length(inInd2)==0
    % ���û��Ѫ���������㣬˵���ö�Ѫ����һ������߽�
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
    % �����������������Ϊ�Ӻ�Ѫ��
    BCType=[BCType 'J'];
    Bifur(i,1:2)=inInd2;
  end
  
  %%% ��������߽�
  outInd1=find(To(i)==From);  % ��ѯ�ö�Ѫ�������Ѫ��
  outInd2=find(To(i)==To);    % ��ѯ��������ö�Ѫ���յ��Ѫ��
  if length(outInd1)==2
    % ����ö�Ѫ����������Ѫ�ܣ�˵���ýڵ���һ���ֲ�ڵ�
    BCType=[BCType 'B'];
    Bifur(i,3:4)=outInd1;
  elseif length(outInd2)==2
    % �����2��Ѫ�ܵ�����յ㣬˵���ýڵ���һ����۽ڵ�
    BCType=[BCType 'C'];
    outInd2(outInd2==i)=[];
    Bifur(i,3:4)=[outInd2 outInd1];
  elseif length(outInd1)==0
    % ������յ㲻���κ�Ѫ�ܵ���㣬˵���ýڵ���һ�������߽�
    % ��̬ģ��
    BCType=[BCType 'p'];
    BCVal(i)=SS_Press(i);
%   else
%     % ����߽�����Ϊ����������������Ϊ�߽磬����Ϊ��
%     BCType=[BCType 'q'];
%     BCVal(i)=-40000/60/1e12;
  else
    % �����������������Ϊ�Ӻ�Ѫ��
    BCType=[BCType 'J'];
    Bifur(i,3:4)=outInd1;
  end
  BCTypeAll=[BCTypeAll;BCType];
end