function [AdapPara,Ev]=AdapLinEqu(NetTypeID,DampFactor,VesNum,SegType,From,To,Diam,Len,Visc,WallTh,Boundary,AdapType,FuncPara,DataArray)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN E5 E02 Visc5 Visc02 DampE DampVisc

% 边界数据
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% 边界血管方向修正 %%%%
% 依据边界条件的符号修正
% 血流为正数->入边界->边界节点在From类中
% 血流为负数->出边界->边界节点在To类中
[From,To]=OriNodeModify(Boundary,From,To);

% 单位调整
% NOTICE: 与其它子程序不同，此处单位调整成mm，参照王若帆的程序
Len=Len*1e-3;   %mm
Diam=Diam*1e-3;   %mm
Visc=Visc*1e-3;   %Pa.s
WallTh=WallTh*1e-3; %mm

% 自适应参数，初值
Qref=0.001;    %流量修正参数，对流模块
PO2ref=94.4;   %氧分压对照参数，对流模块
Tauref=0.095;  %剪切力修正参数，剪切力模块
Jo=7142.8;     %传导信号计算参数，传导模块
Lref=24530;    %血管长度衰减参数，传导模块
t=0.2;
Mo=1000;

% x=[kp,km,kc,ks];

if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  % if 0
  % 如果有自适应数据，那么只跑一次
  Loop1_Visc_Num=1;
  Loop2_Adap_Num=1;
  % 读取FuncPara里的数据作为边界
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo,FuncPara);
else
  % 最大循环次数 TODO(panqing):以收敛判断
  Loop1_Visc_Num=50;
  Loop2_Adap_Num=1000;
  % 无FuncPara
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo);
end

% 记录每次循环得到的参数
DebugVisc=zeros(VesNum,Loop1_Visc_Num);
DebugHd=zeros(VesNum,Loop1_Visc_Num);
DebugP=zeros(VesNum,Loop1_Visc_Num);
DebugFlow=zeros(VesNum,Loop1_Visc_Num);

% 初始Hd计算顺序. Porder, 正序. Norder, 逆序
Porder=1:VesNum;
Norder=VesNum:-1:1;

% Simplex Downhill optimization
options=optimset('tolfun',1e-8,'tolx',1e-4,'MaxIter',10,'Display','iter');
[X,FVAL,EXITFLAG,OUTPUT]= fminsearch(@(x) AdapFunc_WallTh(x,Porder,Norder,Qref,PO2ref,Tauref,Jo,Lref,t,Mo,...
  Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,WallTh,...
  From,To,BHd,BSO2in,BJm,BJc,DataArray),[1.66;0.955;-0.374;3.077;0.0177;0.114;0.609],options);
Tauref=0.095;
Jo=7142.8;
Lref=24530;
kp=X(1);
km=X(2);
kc=X(3);
ks=X(4);
Qref=0.001;
PO2ref=94.4;
AdapPara=[kp,km,kc,ks,PO2ref,Qref,Tauref,Lref,Jo];
Ev=FVAL;
% % 使用优化后的参数进行一次自适应运算
% Ev=AdapFunc(X,Porder,Norder,Qref,PO2ref,Tauref,Jo,Lref,t,Mo,...
%   Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,...
%   From,To,BHd,BSO2in,BJm,BJc,DataArray);