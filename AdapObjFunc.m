function [ObjFuncValue,AdapCoeff,OptDiam,OptWallTh,PO2,Sm,Sc] = AdapObjFunc(AdapCoeff,NetTypeID,AdapType,HdOrder,...
  AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray)
global WITH_WALL WITHOUT_WALL WITH_WALL_Cx WITHOUT_WALL_Cx WITH_WALL_Cx_PULSE

%% 1. �����Ԥ��������
% 1.1 ���Ѫ����������
From=DatMatrix(:,2);
To=DatMatrix(:,3);
Len=DatMatrix(:,4);
Diam=DataArray(:,5);
WallTh=DatMatrix(:,6);
Visc=DatMatrix(:,8);
VesNum=length(From);
% NOTICE: �������ӳ���ͬ���˴���λ������mm�������������ĳ���
Len=Len*1e-3;       % mm
Diam=Diam*1e-3;     % mm
Visc=Visc*1e-3;     % Pa.s
% WallTh=WallTh*1e-3; % mm

% 1.2 ����߽�����
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);
BHd=AdapBoundary(:,1);
BSO2in=AdapBoundary(:,2);
BJm=AdapBoundary(:,3);
BJc=AdapBoundary(:,4);

% 1.3 �������Ӧ�������
dt=AdapPara(1);
M0=AdapPara(2);
Qref=AdapPara(3);
PO2ref=AdapPara(4);
J0=AdapPara(5);

% 1.4 Hd����˳��
Porder=HdOrder(1,:);
Norder=HdOrder(2,:);

% 1.5 ����ԭʼDiam��Visc
OrgDiam=Diam;

%% 2.��������Ӧ��������
% 2.1 ��ʼ����������
Loop1_Visc_Limit=30;
Loop2_Adap_Limit=1000;

% 2.2 ��������Ӧ����
switch AdapType
  case WITH_WALL
    % Ҫ�Ż�������Ӧϵ��
    kc=AdapCoeff(1);
    kmd=AdapCoeff(2);
    kmg=AdapCoeff(3);
    ksd=AdapCoeff(4);
    ksg=AdapCoeff(5);
    kwt=AdapCoeff(6);
    kwo=AdapCoeff(7);
    tauref=AdapCoeff(8);
    J0=AdapCoeff(9);
    Lref=AdapCoeff(10);
    Oref=AdapCoeff(11);
    wref=AdapCoeff(12);
  case WITHOUT_WALL
    % TODO
  case {WITH_WALL_Cx,WITH_WALL_Cx_PULSE}
    % ����Ԥ�������Ӧϵ��
    load AdapCoeff_Wall.mat;
    kc=AdapCoeff_Wall(1);
    kmd=AdapCoeff_Wall(2);
    kmg=AdapCoeff_Wall(3);
    ksd=AdapCoeff_Wall(4);
    ksg=AdapCoeff_Wall(5);
    kwt=AdapCoeff_Wall(6);
    kwo=AdapCoeff_Wall(7);
    tauref=AdapCoeff_Wall(8);
    J0=AdapCoeff_Wall(9);
    Lref=AdapCoeff_Wall(10);
    Oref=AdapCoeff_Wall(11);
    wref=AdapCoeff_Wall(12);
    % Ҫ�Ż���Cxϵ��
    Cbasal=AdapCoeff(1);
    Cpeak=AdapCoeff(2);
    tau_peak=AdapCoeff(3);
    tau_width=AdapCoeff(4);
  case WITHOUT_WALL_Cx
    % TODO
end

%% 3. ��ʼ����Ӧѭ��
% ��ʼ��ѭ��������
Loop2_Cnt=0;
while Loop2_Cnt<Loop2_Adap_Limit % ѭ��2������Ӧѭ������
  Loop2_Cnt=Loop2_Cnt+1;
  Loop1_Cnt=0;
  Loop1_Visc_Num=Loop1_Visc_Limit;
  % ���ݼ�¼����
  DebugVisc=zeros(VesNum,Loop1_Visc_Num);
  FromNew=From;ToNew=To;
  while Loop1_Cnt<Loop1_Visc_Num   % ѭ��1��Visc����
    Loop1_Cnt=Loop1_Cnt+1;
    % ���Է������ģ��
    [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam/1e3,Visc,Len/1e3,FromNew,ToNew);
    % ��������ģ��
    [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(FromNew,ToNew,MeanFlow,DeltaP,Eju);
    % ����Hd����˳��
    [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,1,Eju);
    if Eju==1
      errorFlag=1;
      break;
    end
    % ����Hd
    [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam/1e3,MeanFlowNew,Eju);
    
    %%ճ�Ͷȷ���
    umDiam=Diam.*1e3;   % um
    for i=1:VesNum
      Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
    end
    DebugVisc(:,Loop1_Cnt)=Visc;
    % ѭ������ڶ��κ�(Loop1_Cnt>1)����ʼУ��Visc
    % ����SOR����У�����������ճ�Ͷ��񵴵����
    ViscModType=1;  % 0����������1����Bifurcation������2����All segments����
    Alpha=0.5;
    DebugVisc=modifyViscosity(Porder,FromNew,ToNew,DebugVisc,Loop1_Cnt,ViscModType,Alpha,Eju);
    Visc=DebugVisc(:,Loop1_Cnt);
    
    %%%% �ж�ѭ������ģ�� %%%%
    % Loop1OutType:0-������1-������2-�ﵽ����δ������3-visc����4-���Է��̴���
    [ViscMAE,Loop1OutType,Loop1_Visc_Num]=ViscLoopTerminator(DebugVisc*1e3,Loop1_Visc_Num,Loop1_Cnt,Eju);
  end
  
  % Ѫ������ѧѭ�������󣬼�������Ӧ����
  switch AdapType
    case {WITH_WALL,WITH_WALL_Cx,WITH_WALL_Cx_PULSE}
      % Ref. Pries AR, Reglin B, and Secomb TW. Remodeling of blood vessels:
      % Responses of diameter and wall thickness to hemodynamic and metabolic stimuli.
      % Hypertension 46: 725-731, 2005.
      Dm=Diam+WallTh/2;       % �бڹܾ�(mm)
      Aw=WallTh.*2.*pi.*Dm;   % �бڽ����(mm)
      w=WallTh*1e3;           % �ں�(um)
      tau=32.*Visc.*MeanFlowNew/60*1e-12./(pi*Diam.^3)*1e9*10;  % ������ dyne/cm2;
      %���¼��㹫ʽ�У�DeltaP:mmHg; Diam:mm; WallTh:mm
      O=133.*abs(MeanP).*Diam./(2.*WallTh)*10;                  % ��Ӧ��  dyne/cm2
      if AdapType==WITH_WALL_Cx || AdapType==WITH_WALL_Cx_PULSE
        Cex=Cx40exp(tau,Cbasal,Cpeak,tau_peak,tau_width);
        LrefVec=Lref*Cex.^0.5;
      else
        LrefVec=Lref.*ones(VesNum,1);
      end
      
      % ����������Ӧ�źż���
      % PO2����ģ��
      [PO2,SO2in,SO2mid,SO2out]=PO2Counter(Porder,BSO2in,FromNew,ToNew,Hd,Len,MeanFlowNew,Eju);
      % Sm����ģ��
      [Sm,Jm]=SmCounter(Porder,BJm,FromNew,ToNew,Len,MeanFlowNew,PO2,Qref,PO2ref,dt,M0,Eju);
      % Sc����ģ��
      [Sc,Jc]=ScCounter(Norder,BJc,FromNew,ToNew,Len,Sm,J0,LrefVec,Eju);
      
      ktd=1;
      kog=1;
      ee=1e-9;
      Rtd=0.8;Rod=0.05;Rtg=0.05;Rog=0.8;Rw=0.5;
      Stm=ktd.*log10(tau/tauref+ee)./(1+kwt.*log10(w/wref+ee))+kmd.*(Sm+kc.*Sc)-ksd;
      Som=kog.*log10(O/Oref+ee)./(1+kwo.*log10(w/wref+ee))+kmg.*(Sm+kc.*Sc)-ksg;
      dDm=(Rtd.*Stm+Rod.*Som).*Dm.*dt;
      dAw=Rw.*(Rtg.*Stm+Rog.*Som).*Aw.*dt;
      newDm=Dm+dDm;     % mm
      newAw=Aw+dAw;     % mm2
      % ��������WallTh��Diam
      WallTh=newAw./(2.*pi.*newDm);   %mm
      Diam=newDm-WallTh/2;   %mm
      
      if AdapType==WITH_WALL_Cx_PULSE
        % �õ����µ�WallTh��Diam�����ɶ�̬ģ�ͷ����ļ�
        DatMatrix(:,2)=FromNew;
        DatMatrix(:,3)=ToNew;
        DatMatrix(:,5)=Diam*1e3;
        DatMatrix(:,6)=WallTh*1e3;
        DatMatrix(:,8)=Visc;
        SIM_PREP_FUNC(NetTypeID,Boundary,DatMatrix);
        cd('PulseAdapDIR');
        system('run.bat');
        % ���������з���
        cd('..');
      end
    case WITHOUT_WALL
      % TODO
    case WITHOUT_WALL_Cx
      % TODO
  end
  
  %% ��Ҫ������¼
  %   DebugStau(:,Loop2_Cnt)=Stau;
  %   DebugSp(:,Loop2_Cnt)=Sp;
  %   DebugSm(:,Loop2_Cnt)=Sm;
  %   DebugSc(:,Loop2_Cnt)=Sc;
  DebugDiam(:,Loop2_Cnt)=Diam*1e3;
  DebugFlow(:,Loop2_Cnt)=MeanFlowNew;
  DebugP(:,Loop2_Cnt)=MeanP;
  DebugPO2(:,Loop2_Cnt)=PO2;
  %   DebugTau(:,Loop2_Cnt)=Tau;
  DebugSO2in(:,Loop2_Cnt)=SO2in;
  DebugSO2out(:,Loop2_Cnt)=SO2out;
  DebugTHd(:,Loop2_Cnt)=Hd;
  
  %% ����Ӧѭ������ģ�� %%%%
  switch AdapType
    case {WITH_WALL,WITH_WALL_Cx,WITH_WALL_Cx_PULSE}
      if mean(abs(dDm))<1e-5*dt && mean(abs(dAw))<1e-5*dt
        % �ﵽ���ȣ���������Ӧѭ��
        AdapLoopOutType=1;
        break;
      elseif isnan(mean(abs(dDm))) || isnan(mean(abs(dAw)))
        AdapLoopOutType=4;
        break;
      elseif Loop1OutType==2
        break;
      end
      % �����ʾ
      fprintf('AdapIt:=%d DiamMAE=%5.4e dDm:%5.4e dAw=%5.4e ViscIt:%3d.\n',Loop2_Cnt,1,mean(abs(dDm)),mean(abs(dAw)),Loop1_Cnt);
    case {WITHOUT_WALL,WITHOUT_WALL_Cx}
      % AdapLoopOutType:0-������1-������2-�ﵽ����δ������3-visc����4-���Է��̴���
      % AccuracyType�����ȿ��ƣ���0-���ܾ�������1-���ܾ������2-ƽ���ܾ�������3-ƽ���ܾ�������
      AccuracyType=3;
      [DiamMAE,AdapLoopOutType,Loop2_Adap_Num]=AdapLoopTerminator(DebugDiam,Loop2_Adap_Num,Loop2_Cnt,Eju,AccuracyType);
      % �쳣���������������
      if AdapLoopOutType>2
        if Loop2_Cnt==1
          Diam=OrgDiam;
        else  % ��������ʱ�����ǰһ��ѭ�����ݣ��Ա�����������ļ���
          % Ϊ�˼���Ev��Ed
          Diam=DebugDiam(:,Loop2_Cnt-1)./1e3;
        end
      end
      % �����ʾ
      fprintf('AdapIt:=%d DiamMAE=%5.4e dDm:%5.4e dAw=%5.4e ViscIt:%3d.\n',Loop2_Cnt,1,mean(abs(dDm)),mean(abs(dAw)),Loop1_Cnt,ViscMAE);
  end
end

%% 4. ����Ӧ�����Ȳ����Ż����۲��� %%%%
% Ed���ܾ��仯���죻Ev�����ٱ仯����
umDiam=Diam*1e3;
% ���ټ���
Vel=4.*1e3.*MeanFlowNew./(60.*pi.*umDiam.*umDiam);   %Ѫ������mm/s
if isempty(DataArray)
  Ed=sqrt(mean((OrgDiam-Diam).^2./(OrgDiam.^2)));
else
  Ev=sqrt(mean(4*(abs(DataArray(:,8))-Vel).^2./(abs(DataArray(:,8))+Vel).^2));
  Ed=sqrt(mean(4*(DataArray(:,5)-umDiam).^2./(DataArray(:,5).^2)));
end
% ���ڷ�����������Ev�ͷ�
if Loop1OutType~=1 || Loop2_Cnt == 1000 || AdapLoopOutType~=1
  Ev=Ev*10;
  Ed=Ed*10;
end
% ��¼�Ż�Ŀ�꺯��ֵ����Ӧ���Ż����
ObjFuncValue=Ed;
OptDiam=Diam;
OptWallTh=WallTh;
fprintf('ObjFuncValue=%5.4e\n',ObjFuncValue);
