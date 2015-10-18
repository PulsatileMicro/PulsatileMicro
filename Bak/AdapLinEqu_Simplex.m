function [AdapPara,Ev]=AdapLinEqu(NetTypeID,DampFactor,VesNum,SegType,From,To,Diam,Len,Visc,WallTh,Boundary,AdapType,FuncPara,DataArray)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN E5 E02 Visc5 Visc02 DampE DampVisc

% �߽�����
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% �߽�Ѫ�ܷ������� %%%%
% ���ݱ߽������ķ�������
% Ѫ��Ϊ����->��߽�->�߽�ڵ���From����
% Ѫ��Ϊ����->���߽�->�߽�ڵ���To����
[From,To]=OriNodeModify(Boundary,From,To);

% ��λ����
% NOTICE: �������ӳ���ͬ���˴���λ������mm�������������ĳ���
Len=Len*1e-3;   %mm
Diam=Diam*1e-3;   %mm
Visc=Visc*1e-3;   %Pa.s
WallTh=WallTh*1e-3; %mm

% ����Ӧ��������ֵ
Qref=0.001;    %������������������ģ��
PO2ref=94.4;   %����ѹ���ղ���������ģ��
Tauref=0.095;  %����������������������ģ��
Jo=7142.8;     %�����źż������������ģ��
Lref=24530;    %Ѫ�ܳ���˥������������ģ��
t=0.2;
Mo=1000;

% x=[kp,km,kc,ks];

if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  % if 0
  % ���������Ӧ���ݣ���ôֻ��һ��
  Loop1_Visc_Num=1;
  Loop2_Adap_Num=1;
  % ��ȡFuncPara���������Ϊ�߽�
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo,FuncPara);
else
  % ���ѭ������ TODO(panqing):�������ж�
  Loop1_Visc_Num=50;
  Loop2_Adap_Num=1000;
  % ��FuncPara
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo);
end

% ��¼ÿ��ѭ���õ��Ĳ���
DebugVisc=zeros(VesNum,Loop1_Visc_Num);
DebugHd=zeros(VesNum,Loop1_Visc_Num);
DebugP=zeros(VesNum,Loop1_Visc_Num);
DebugFlow=zeros(VesNum,Loop1_Visc_Num);

% ��ʼHd����˳��. Porder, ����. Norder, ����
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
% % ʹ���Ż���Ĳ�������һ������Ӧ����
% Ev=AdapFunc(X,Porder,Norder,Qref,PO2ref,Tauref,Jo,Lref,t,Mo,...
%   Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,...
%   From,To,BHd,BSO2in,BJm,BJc,DataArray);