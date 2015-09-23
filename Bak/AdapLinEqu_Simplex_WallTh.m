function [AdapPara,Ev,WallTh,Diam,Visc,Sm,Sc]=AdapLinEqu(NetTypeID,DampFactor,VesNum,SegType,From,To,Diam,Len,Visc,WallTh,Boundary,AdapType,FuncPara,DataArray)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN E5 E02 Visc5 Visc02 DampE DampVisc

% 边界血管方向修正 %%%%
% 依据边界条件的符号修正
% 血流为正数->入边界->边界节点在From类中
% 血流为负数->出边界->边界节点在To类中
[From,To]=OriNodeModify(Boundary,From,To);

% 单位调整
% NOTICE: 与其它子程序不同，此处单位调整成mm，参照王若帆的程序
Len=Len*1e-3;       %mm
Diam=Diam*1e-3;     %mm
Visc=Visc*1e-3;     %Pa.s
WallTh=WallTh*1e-3; %mm

% 自适应参数，初值
Qref=0.001;    % 流量修正参数，对流模块
PO2ref=94.4;   % 氧分压对照参数，对流模块
Jo=6618;       % 传导信号计算参数，传导模块
t=0.1;
Mo=1000;

if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
% if 0
  % 如果有自适应数据，那么只跑一次
  Loop1_Visc_Num=30;
  Loop2_Adap_Num=1000;
  % 读取FuncPara里的数据作为边界
  WallTh=4e-3*ones(VesNum,1);
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo,FuncPara);
else
  % 最大循环次数 TODO(panqing):以收敛判断
  Loop1_Visc_Num=30;
  Loop2_Adap_Num=1000;
  % 无FuncPara
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo);
end

% 初始Hd计算顺序. Porder, 正序. Norder, 逆序
Porder=1:VesNum;
Norder=VesNum:-1:1;

% WallAdapPara=[1.66,0.955,-0.374,3.077,0.0177,0.114,0.609,0.5598,6618,14292,32050,0.804];
load AdapWallPara.mat;
if AdapType==0
  load AdapWallPara.mat;
  AdapPara=WallAdapPara;
  % 使用优化后的参数进行一次自适应运算
  [Ev,WallTh,Diam,Visc,Sm,Sc]=AdapFunc_WallTh(AdapPara,Porder,Norder,Qref,PO2ref,t,Mo,...
    Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,WallTh,...
    From,To,BHd,BSO2in,BJm,BJc,DataArray,AdapType,[]);
elseif AdapType==1
  AdapPara=WallAdapPara;
  % Simplex Downhill optimization
  options=optimset('tolfun',1e-3,'tolx',1e-2,'MaxIter',50,'Display','iter');
  [X,FVAL,EXITFLAG,OUTPUT]= fminsearch(@(x) AdapFunc_WallTh(x,Porder,Norder,Qref,PO2ref,t,Mo,...
    Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,WallTh,...
    From,To,BHd,BSO2in,BJm,BJc,DataArray,AdapType,[]),AdapPara,options);
  AdapPara=X;
  % 使用优化后的参数进行一次自适应运算
  [Ev,WallTh,Diam,Visc,Sm,Sc]=AdapFunc_WallTh(AdapPara,Porder,Norder,Qref,PO2ref,t,Mo,...
    Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,WallTh,...
    From,To,BHd,BSO2in,BJm,BJc,DataArray,AdapType,[]);
elseif AdapType==2
  CxAdapPara=[0.5,3.31,41.1,55.6];
  % Simplex Downhill optimization
  options=optimset('tolfun',1e-3,'tolx',1e-2,'MaxIter',50,'Display','iter');
  [X,FVAL,EXITFLAG,OUTPUT]= fminsearch(@(x) AdapFunc_WallTh(x,Porder,Norder,Qref,PO2ref,t,Mo,...
    Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,WallTh,...
    From,To,BHd,BSO2in,BJm,BJc,DataArray,AdapType,WallAdapPara),CxAdapPara,options);
    AdapPara=X;
    % 使用优化后的参数进行一次自适应运算
  [Ev,WallTh,Diam,Visc,Sm,Sc]=AdapFunc_WallTh(AdapPara,Porder,Norder,Qref,PO2ref,t,Mo,...
    Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,WallTh,...
    From,To,BHd,BSO2in,BJm,BJc,DataArray,AdapType,WallAdapPara);
end

WallTh=WallTh*1e3;
Diam=Diam*1e3;