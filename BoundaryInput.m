function [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo,FuncPara)
if nargin<5
  % 只有血流动力学
  DataType=2;
elseif nargin<6
  % 没有FuncPara
  DataType=1;
else
  % 有FuncPara
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
        % 从FuncPara中读取数据
        BJm(FIndex)=JmBoundary(FuncPara(FIndex,8),FuncPara(FIndex,14),Qref);
      end
      if ~isempty(TIndex)
        % 从FuncPara中读取数据
        BJc(TIndex)=JcBoundary(FuncPara(TIndex,16),Jo);
      end
    end           % 边界条件中的Hd赋值
  case 1
    % 没有FuncPara
    for i=1:length(Boundary(:,1))
      FIndex=find(From==Boundary(i,1));
      TIndex=To==Boundary(i,1);
      BHd(FIndex)=Boundary(i,4);
      BSO2in(FIndex)=0.94;
      BJm(FIndex)=0.01;    % 该如何设置？
      BJc(TIndex)=0.01;
    end           % 边界条件中的Hd赋值
  case 2
    % 只有血流动力学部分
    for i=1:length(Boundary(:,1))
      FIndex=find(From==Boundary(i,1));
      TIndex= To==Boundary(i,1);
      BHd(FIndex)=Boundary(i,4);
    end           %边界条件中的Hd赋值
end