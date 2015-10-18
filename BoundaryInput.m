function [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo,FuncPara)
if nargin<5
  % ֻ��Ѫ������ѧ
  DataType=2;
elseif nargin<6
  % û��FuncPara
  DataType=1;
else
  % ��FuncPara
  DataType=0;
end
VesNum=length(From);
BHd=zeros(VesNum,1);
BSO2in=zeros(VesNum,1);
BJm=zeros(VesNum,1);
BJc=zeros(VesNum,1);
switch DataType
  case 0
    for i=1:length(Boundary(:,1))
      FIndex=find(From==Boundary(i,1));
      TIndex=find(To==Boundary(i,1));
      if ~isempty(FIndex)
        BHd(FIndex)=Boundary(i,4);
        BSO2in(FIndex)=FuncPara(FIndex,22);
        % ��FuncPara�ж�ȡ����
        BJm(FIndex)=JmBoundary(FuncPara(FIndex,8),FuncPara(FIndex,14),Qref);
      end
      if ~isempty(TIndex)
        % ��FuncPara�ж�ȡ����
        BJc(TIndex)=JcBoundary(FuncPara(TIndex,16),Jo);
      end
    end           % �߽������е�Hd��ֵ
  case 1
    % û��FuncPara
    for i=1:length(Boundary(:,1))
      FIndex=find(From==Boundary(i,1));
      TIndex=To==Boundary(i,1);
      BHd(FIndex)=Boundary(i,4);
      BSO2in(FIndex)=0.94;
      BJm(FIndex)=0.01;    % ��������ã�
      BJc(TIndex)=0.01;
    end           % �߽������е�Hd��ֵ
  case 2
    % ֻ��Ѫ������ѧ����
    for i=1:length(Boundary(:,1))
      FIndex=find(From==Boundary(i,1));
      TIndex= To==Boundary(i,1);
      BHd(FIndex)=Boundary(i,4);
    end           %�߽������е�Hd��ֵ
end