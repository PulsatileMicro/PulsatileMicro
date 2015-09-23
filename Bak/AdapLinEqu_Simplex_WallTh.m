function [AdapPara,Ev,WallTh,Diam,Visc,Sm,Sc]=AdapLinEqu(NetTypeID,DampFactor,VesNum,SegType,From,To,Diam,Len,Visc,WallTh,Boundary,AdapType,FuncPara,DataArray)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN E5 E02 Visc5 Visc02 DampE DampVisc

% �߽�Ѫ�ܷ������� %%%%
% ���ݱ߽������ķ�������
% Ѫ��Ϊ����->��߽�->�߽�ڵ���From����
% Ѫ��Ϊ����->���߽�->�߽�ڵ���To����
[From,To]=OriNodeModify(Boundary,From,To);

% ��λ����
% NOTICE: �������ӳ���ͬ���˴���λ������mm�������������ĳ���
Len=Len*1e-3;       %mm
Diam=Diam*1e-3;     %mm
Visc=Visc*1e-3;     %Pa.s
WallTh=WallTh*1e-3; %mm

% ����Ӧ��������ֵ
Qref=0.001;    % ������������������ģ��
PO2ref=94.4;   % ����ѹ���ղ���������ģ��
Jo=6618;       % �����źż������������ģ��
t=0.1;
Mo=1000;

if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
% if 0
  % ���������Ӧ���ݣ���ôֻ��һ��
  Loop1_Visc_Num=30;
  Loop2_Adap_Num=1000;
  % ��ȡFuncPara���������Ϊ�߽�
  WallTh=4e-3*ones(VesNum,1);
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo,FuncPara);
else
  % ���ѭ������ TODO(panqing):�������ж�
  Loop1_Visc_Num=30;
  Loop2_Adap_Num=1000;
  % ��FuncPara
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo);
end

% ��ʼHd����˳��. Porder, ����. Norder, ����
Porder=1:VesNum;
Norder=VesNum:-1:1;

% WallAdapPara=[1.66,0.955,-0.374,3.077,0.0177,0.114,0.609,0.5598,6618,14292,32050,0.804];
load AdapWallPara.mat;
if AdapType==0
  load AdapWallPara.mat;
  AdapPara=WallAdapPara;
  % ʹ���Ż���Ĳ�������һ������Ӧ����
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
  % ʹ���Ż���Ĳ�������һ������Ӧ����
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
    % ʹ���Ż���Ĳ�������һ������Ӧ����
  [Ev,WallTh,Diam,Visc,Sm,Sc]=AdapFunc_WallTh(AdapPara,Porder,Norder,Qref,PO2ref,t,Mo,...
    Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,WallTh,...
    From,To,BHd,BSO2in,BJm,BJc,DataArray,AdapType,WallAdapPara);
end

WallTh=WallTh*1e3;
Diam=Diam*1e3;