function [y,WallTh,Diam,Visc,Sm,Sc]= AdapFunc_WallTh(X,Porder,Norder,Qref,PO2ref,t,Mo,...
  Loop2_Adap_Num,Loop1_Visc_Num_Init,VesNum,Boundary,Diam,Len,Visc,WallTh,...
  From,To,BHd,BSO2in,BJm,BJc,DataArray,AdapType,WallAdapPara)

% �߽�����
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% ����ԭʼDiam��Visc
OrgDiam=Diam;
OrgVisc=Visc;

Loop2_Cnt=0;
% Opt for Wall thickness
if AdapType~=2
  WallAdapPara=X;
else
  % Opt for Cx
  Cbasal=X(1);
  Cpeak=X(2);
  tau_peak=X(3);
  tau_width=X(4);
end
kc=WallAdapPara(1);
kmd=WallAdapPara(2);
kmg=WallAdapPara(3);
ksd=WallAdapPara(4);
ksg=WallAdapPara(5);
kwt=WallAdapPara(6);
kwo=WallAdapPara(7);
tauref=WallAdapPara(8);
Jo=WallAdapPara(9);
Lref=WallAdapPara(10);
Oref=WallAdapPara(11);
wref=WallAdapPara(12);
LrefVec=Lref*ones(VesNum,1);

while Loop2_Cnt<Loop2_Adap_Num % ѭ��2������Ӧ����
  Loop2_Cnt=Loop2_Cnt+1;
  Loop1_Cnt=0;
  Loop1_Visc_Num=Loop1_Visc_Num_Init;
  % ���ݼ�¼����
  DebugVisc=zeros(VesNum,Loop1_Visc_Num);
  FromNew=From;ToNew=To;
  while Loop1_Cnt<Loop1_Visc_Num   % ѭ��1��Visc����  
    Loop1_Cnt=Loop1_Cnt+1;
    % ���Է������ģ��
    [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam/1e3,Visc,Len/1e3,FromNew,ToNew);
    % ��������ģ��
    [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(FromNew,ToNew,MeanFlow,DeltaP,Eju);
    %   FromNew=From;ToNew=To;MeanFlowNew=MeanFlow;DeltaPNew=DeltaP;
    % ����Hd����˳��
    [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,2,Eju);
    if Eju==1
      errorFlag=1;
      break;
    end
    % ����Hd
    [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam/1e3,MeanFlowNew,Eju);
    
    % ����Ӧ���������治��Ҫճ�Ͷȵ���
    %       if NetTypeID~=Net_546_ID && NetTypeID~=Net_122_ID
    if 1
      %%ճ�Ͷȷ���
      umDiam=Diam.*1e3;   %um
      for i=1:VesNum
        Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
      end
      DebugVisc(:,Loop1_Cnt)=Visc;
      % ѭ������ڶ��κ�(Loop1_Cnt>1)����ʼУ��Visc
      % ����SOR����У�����������ճ�Ͷ��񵴵����
      ViscModType=1;  %0����������1����Bifurcation������2����All segments����
      Alpha=0.5;
      DebugVisc=modifyViscosity(Porder,FromNew,ToNew,DebugVisc,Loop1_Cnt,ViscModType,Alpha,Eju);
      Visc=DebugVisc(:,Loop1_Cnt);
    end
    
    %% ѭ������ģ�� %%%%
    %MAE: max absolute error
    %LoopOutType:0-������1-������2-�ﵽ����δ������3-visc����4-���Է��̴���
    [ViscMAE,Loop1OutType,Loop1_Visc_Num]=ViscLoopTerminator(DebugVisc*1e3,Loop1_Visc_Num,Loop1_Cnt,Eju);
  end
  
  % ����Ӧ����
  % Hypertension���׷���
  Dm=Diam+WallTh/2;    %�бڹܾ�(mm)
  Aw=WallTh.*2.*pi.*Dm; %�бڽ����(mm)
  w=WallTh*1e3;   %�ں�(um)
  tau=32.*Visc.*MeanFlowNew/60*1e-12./(pi*Diam.^3)*1e9*10;  % ������ dyn/cm2;
  %O=DeltaP*Diam/(2*WallTh)
  %���¼��㹫ʽ�У�DeltaP:mmHg; Diam:mm; WallTh:mm
  O=133.*abs(MeanP).*Diam./(2.*WallTh)*10;    % ��Ӧ��  dyn/cm2
  if AdapType==2
    Cex=Cx40exp(tau,Cbasal,Cpeak,tau_peak,tau_width);
    LrefVec=Lref*Cex.^0.5;
  end
  % PO2����ģ��
  [PO2,SO2in,SO2mid,SO2out]=PO2Counter(Porder,BSO2in,FromNew,ToNew,Hd,Len,MeanFlowNew,Eju);
  % Sm����ģ��
  [Sm,Jm]=SmCounter(Porder,BJm,FromNew,ToNew,Len,MeanFlowNew,PO2,Qref,PO2ref,t,Mo,Eju);
  % Sc����ģ��
  [Sc,Jc]=ScCounter(Norder,BJc,FromNew,ToNew,Len,Sm,Jo,LrefVec,Eju);
%   [Sc,Jc]=ScCounter(Norder,BJc,FromNew,ToNew,Len,Sm,Jo,Lceff,Eju);
  
  % ���� ref.Pries 2005
  ktd=1;
  kog=1;
  ee=1e-9;
  %�����㣺Rtd^2+Rod^2=1, Rtg^2+Rog^2=1
  Rtd=0.8;Rod=0.0;Rtg=0.0;Rog=0.8;Rw=0.5;
  Stm=ktd.*log10(tau/tauref+ee)./(1+kwt.*log10(w/wref+ee))+kmd.*(Sm+kc.*Sc)-ksd;
  Som=kog.*log10(O/Oref+ee)./(1+kwo.*log10(w/wref+ee))+kmg.*(Sm+kc.*Sc)-ksg;
  dDm=(Rtd.*Stm+Rod.*Som).*Dm.*t;
  dAw=Rw.*(Rtg.*Stm+Rog.*Som).*Aw.*t;
  newDm=Dm+dDm;     % mm
  newAw=Aw+dAw;     % mm2
  % ��������WallTh��Diam
  WallTh=newAw./(2.*pi.*newDm);   %mm
  Diam=newDm-WallTh/2;   %mm
  
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
  
  %% ѭ������ģ�� %%%%
  %MAE: max absolute error
  %AdapLoopOutType:0-������1-������2-�ﵽ����δ������3-visc����4-���Է��̴���
  %AccuracyType�����ȿ��ƣ���0-���ܾ�������1-���ܾ������2-ƽ���ܾ�������3-ƽ���ܾ�������
%   AccuracyType=3;
%   [DiamMAE,AdapLoopOutType,Loop2_Adap_Num]=AdapLoopTerminator(DebugDiam,Loop2_Adap_Num,Loop2_Cnt,Eju,AccuracyType);
%   % �쳣���������������
%   if AdapLoopOutType>2
%     if Loop2_Cnt==1
%       Diam=OrgDiam;
%     else  % ��������ʱ�����ǰһ��ѭ�����ݣ��Ա�����������ļ���
%       % Ϊ�˼���Ev��Ed
%       Diam=DebugDiam(:,Loop2_Cnt-1)./1e3;
%     end
%   end    
  
  %% �����ʾ %%%%
  fprintf('Adap Iter:=%d AdapOutType=%d DiamMAE=%5.4e dDm:%5.4e dAw=%5.4e ViscIter:%3d ViscMAE:%5.4e.\n',Loop2_Cnt,1,1,mean(abs(dDm)),mean(abs(dAw)),Loop1_Cnt,ViscMAE);
  if Loop1OutType~=1
    break;
  elseif mean(abs(dDm))<1e-5*t && mean(abs(dAw))<1e-5*t
    break;
  end
end
% ����Ӧ��������¼������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����Ӧ�����Ȳ����Ż����۲��� %%%%
% Ed���ܾ��仯���죻Ev�����ٱ仯����
umDiam=Diam*1e3;
% ���ټ���
Vel=4.*1e3.*MeanFlowNew./(60.*pi.*umDiam.*umDiam);   %Ѫ������mm/s
% Ref. Bettina 2009 (��������
if isempty(DataArray)
  Ed=sqrt(mean((OrgDiam-Diam).^2./(OrgDiam.^2)));
else
  Ev=sqrt(mean(4*(abs(DataArray(:,8))-Vel).^2./(abs(DataArray(:,8))+Vel).^2));
  Ed=sqrt(mean(4*(DataArray(:,5)-Diam*1e3).^2./(DataArray(:,5).^2)));
end
% ���ڷ�����������Ev�ͷ�
if Loop1OutType~=1
  Ev=Ev*10;
  Ed=Ed*10;
end
% if AdapLoopOutType>1
%   if AdapLoopOutType==2
%     Ev=Ev*1.5;
%     Ed=Ed*1.5;
%   else
%     Ev=Ev*2;
%     Ed=Ed*2;
%   end
% end

y=Ev;
fprintf('Ev=%5.4e\n',y);
end
