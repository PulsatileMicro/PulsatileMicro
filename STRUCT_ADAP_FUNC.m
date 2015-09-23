function [NewAdapCoeff,ObjValue,OptDiam,OptWallTh,Sm,Sc] = STRUCT_ADAP_FUNC(NetTypeID,DataArray,Boundary,FuncPara,DatMatrix)
global WITH_WALL WITHOUT_WALL WITH_WALL_Cx WITHOUT_WALL_Cx WITH_WALL_Cx_PULSE
global PSO DOWNHILL GLOBAL_SEARCH
global NOT_OPT_PARA OPT_PARA
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID

%% 1. 设置自适应方法、优化方法、优化类型
AdapType=WITH_WALL;
%%OptMethod=DOWNHILL;
OptMethod=PSO;
OptType=OPT_PARA;

%% 2. 数据预处理
% 2.1 数据解包
From=DatMatrix(:,2);
To=DatMatrix(:,3);
VesNum=length(From);

% 2.2 单位调整
% NOTICE: 与其它子程序不同，此处单位调整成mm，参照王若帆的程序
% Len=Len*1e-3;       % mm
% Diam=Diam*1e-3;     % mm
% Visc=Visc*1e-3;     % Pa.s
% WallTh=WallTh*1e-3; % mm

% 2.3 边界方向修正
% 依据边界条件的符号修正
% 血流为正数->入边界->边界节点在From类中
% 血流为负数->出边界->边界节点在To类中
[From,To]=OriNodeModify(Boundary,From,To);

%% 3. 根据所设参数开始自适应
% 3.1 自适应参数设置
dt=0.1;       % 自适应步长
M0=1000;
Qref=0.001;   % 流量修正参数，对流模块，无需优化
PO2ref=94.4;   % 氧分压对照参数，对流模块
J0=6618;       % 传导信号计算参数，传导模块
AdapPara=[dt,M0,Qref,PO2ref,J0];

switch AdapType
  case WITH_WALL    % 12个参数
    DatMatrix(:,6)=4e-3*ones(VesNum,1);
    load AdapCoeff_Wall.mat;
%     load AdapCoeff_Wall_Pries.mat;
    AdapCoeff=AdapCoeff_Wall;
  case WITHOUT_WALL
    % TODO
  case {WITH_WALL_Cx,WITH_WALL_Cx_PULSE}
    DatMatrix(:,6)=4e-3*ones(VesNum,1);
    load AdapCoeff_Wall.mat;
    load AdapCoeff_Wall_Cx.mat;
    AdapCoeff=AdapCoeff_Wall_Cx;
  case WITHOUT_WALL_Cx
    % TODO
end

if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  % 读取FuncPara里的数据作为边界
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,J0,FuncPara);
else
  % 无FuncPara
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,J0);
end
AdapBoundary=[BHd,BSO2in,BJm,BJc];

% 初始Hd计算顺序. Porder, 正序. Norder, 逆序
Porder=1:VesNum;
Norder=VesNum:-1:1;
HdOrder=[Porder;Norder];

if OptType==OPT_PARA
  % 选择优化算法
  switch OptMethod
    case PSO
      % TODO: 仿真出错的情况需要处理
      SwarmSize=5;
      InitSwarm=zeros(SwarmSize,length(AdapCoeff));
      for i=1:SwarmSize
        InitSwarm(i,:)=AdapCoeff.*(rand(1,length(AdapCoeff))+0.5);
      end
      options = optimoptions('particleswarm','SwarmSize',SwarmSize,'InitialSwarm',InitSwarm,'MaxIter',1000);
      %  'OutputFcns',@OptOutFunc,'PlotFcns',@OptPlotFunc,
      lb=[];
      ub=[];
      [X,FVAL,EXITFLAG,OUTPUT]= particleswarm(@(x) AdapObjFunc(x,NetTypeID,AdapType,HdOrder,...
        AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray),length(AdapCoeff),lb,ub,options);
      NewAdapCoeff=X;
    case DOWNHILL
      % Simplex Downhill方法
      options=optimset('TolFun',1e-2,'TolX',1e-2,'MaxIter',100,'Display','iter');
      [X,FVAL,EXITFLAG,OUTPUT]= fminsearch(@(x) AdapObjFunc(x,NetTypeID,AdapType,HdOrder,...
        AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray),AdapCoeff,options);
      NewAdapCoeff=X;
%   EXITFLAG 
%   1   Relative change in the objective value over the last options.StallIterLimit iterations is less than options.TolFun.
%   0   Number of iterations exceeded options.MaxIter.
%   -1  Iterations stopped by output function or plot function.
%   -2  Bounds are inconsistent: for some i, lb(i) > ub(i).
%   -3  Best objective function value is at or below options.ObjectiveLimit.
%   -4  Best objective function value did not change within options.StallTimeLimit seconds.
%   -5  Run time exceeded options.MaxTime seconds.
      formatstring = 'particleswarm reached the value %f using %d function evaluations, after %d iterations with ExitFlag %d and algorithm stop reason is %s.\n';
      fprintf(formatstring, FVAL, OUTPUT.funccount, OUTPUT.iterations, EXITFLAG, OUTPUT.message);
    case GLOBAL_SEARCH
      % TODO
  end
  % 使用更新后的系数跑一次
  [ObjValue,AdapCoeff,OptDiam,OptWallTh,Sm,Sc]=AdapObjFunc(NewAdapCoeff,NetTypeID,AdapType,HdOrder,...
    AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray);
  NewAdapCoeff=AdapCoeff;
else
  [ObjValue,AdapCoeff,OptDiam,OptWallTh,Sm,Sc]=AdapObjFunc(AdapCoeff,NetTypeID,AdapType,HdOrder,...
    AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray);
  NewAdapCoeff=AdapCoeff;
end