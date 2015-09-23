function y = AdapFunc(X,Porder,Norder,Qref,PO2ref,Tauref,Jo,Lref,t,Mo,...
  Loop2_Adap_Num,Loop1_Visc_Num,VesNum,Boundary,Diam,Len,Visc,WallTh,...
  From,To,BHd,BSO2in,BJm,BJc,DataArray)

% �߽�����
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

OrgDiam=Diam;
OrgVisc=Visc;

Loop2_Cnt=0;
kp=X(1);
km=X(2);
kc=X(3);
ks=X(4);
while Loop2_Cnt<Loop2_Adap_Num % ѭ��2������Ӧ���ؼ���
  Loop2_Cnt=Loop2_Cnt+1;
  Loop1_Cnt=0;
  Loop1_Visc_Num=50;
  
  % ���ݼ�¼����
  DebugVisc=zeros(VesNum,Loop1_Visc_Num);
  
  while Loop1_Cnt<Loop1_Visc_Num   % ѭ��1��Visc����
    Loop1_Cnt=Loop1_Cnt+1;
    % ���Է������ģ��
    [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam/1e3,Visc,Len/1e3,From,To);
    % ��������ģ��
    [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju);
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
    
    % ��¼ÿ�ε����ķ�����
    %       DebugVisc(:,Loop1_Cnt)=Visc;
    %       DebugHd(:,Loop1_Cnt)=Hd;
    %       DebugP(:,Loop1_Cnt)=MeanP;
    %       DebugFlow(:,Loop1_Cnt)=MeanFlowNew;
    
    %% ѭ������ģ�� %%%%
    %MAE: max absolute error
    %LoopOutType:0-������1-������2-�ﵽ����δ������3-visc����4-���Է��̴���
    [ViscMAE,LoopOutType,Loop1_Visc_Num]=ViscLoopTerminator(DebugVisc*1e3,Loop1_Visc_Num,Loop1_Cnt,Eju);
  end
  
  %PO2����ģ��
  [PO2,SO2in,SO2mid,SO2out]=PO2Counter(Porder,BSO2in,FromNew,ToNew,Hd,Len,MeanFlowNew,Eju);
  %Sm����ģ��
  [Sm,Jm]=SmCounter(Porder,BJm,FromNew,ToNew,Len,MeanFlowNew,PO2,Qref,PO2ref,t,Mo,Eju);
  %Sc����ģ��
  [Sc,Jc]=ScCounter(Norder,BJc,FromNew,ToNew,Len,Sm,Jo,Lref,Eju);
  %Stau��Sp����ģ��
  % TODO: ��Ҫ�޸�ΪHypertension�����еķ���
  % [Stau,Sp,Tau]=SCounter(DeltaP,Diam,Len,MeanP,Tauref,Eju);
  %�ܾ�����ģ��
  % Stot=Stau+kp.*Sp+km.*(Sm+kc.*Sc)-ks;
  % Diam=Diam+Stot.*Diam.*t;
  
  % Hypertension���׷���
  Dm=Diam+WallTh/2;    %�бڹܾ�(mm)
  Aw=WallTh.*2.*pi.*Dm; %�бڽ����(mm)
  w=WallTh*1e3;   %�ں�(um)
  tau=32*10^4/60.*Visc.*MeanFlowNew./(umDiam.^3);  %������ dyn/cm2;
  %O=DeltaP*Diam/(2*WallTh)
  %���¼��㹫ʽ�У�DeltaP:mmHg; Diam:mm; WallTh:mm
  O=1333.2.*abs(DeltaP).*Diam./(2.*WallTh);    %��Ӧ��  dyn/cm2
  %���� ref.Pries 2005
  ktd=1;
  kog=1;
  tauref=0.5598; %dyn/cm2
  Oref=32050;  %dyn/cm2
  wref=0.804;  %um
  kc=1.66;
  kmd=0.955;
  kmg=-0.374;
  ksd=3.077;
  ksg=0.0177;
  kwt=0.114;
  kwo=0.609;
  ee=0.1;
  %Rtd^2+Rod^2=1
  %Rtg^2+Rog^2=1
  Rtd=0.6;
  Rod=0.8;
  Rtg=0.6;
  Rog=0.8;
  Rw=0.2;
  Stm=ktd.*log10(tau/tauref+ee)./(1+kwt.*log10(w/wref+ee))+kmd.*(Sm+kc.*Sc)-ksd;
  Som=kog.*log10(O/Oref+ee)./(1+kwo.*log10(w/wref+ee))+kmg.*(Sm+kc.*Sc)-ksg;
  dDm=(Rtd.*Stm+Rod.*Som).*Dm.*t;
  dAw=Rw.*(Rtg.*Stm+Rog.*Som).*Aw.*t;
  newDm=Dm+dDm;     %mm
  newAw=Aw+dAw;     %mm2
  %��������
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
  %AccuracyType�����ȿ��ƣ���0-��������1-��������2-ƽ���������
  AccuracyType=2;
  [DiamMAE,AdapLoopOutType,Loop1_Visc_Num]=AdapLoopTerminator(DebugDiam,Loop1_Visc_Num,Loop2_Cnt,Eju,AccuracyType);
  %�쳣���������������
  if AdapLoopOutType>2
    if Loop2_Cnt==1
      Diam=OrgDiam;
    else  %��������ʱ�����ǰһ��ѭ�����ݣ��Ա�����������ļ���
      %Ϊ�˼���Ev��Ed
      Diam=DebugDiam(:,Loop2_Cnt-1)./1e3;
      %Ϊ�˼����Ȩ��·����������ɢ
      MeanFlow=abs(DebugFlow(:,Loop2_Cnt-1));
      % Ϊ�˼���������
      SO2in=DebugSO2in(:,Loop2_Cnt-1);
      SO2out=DebugSO2out(:,Loop2_Cnt-1);
      Hd=DebugTHd(:,Loop2_Cnt-1);
      % Ϊ�˼������ëϸѪ��ѹ
      MeanP=DebugP(:,Loop2_Cnt-1);
    end
  end
  
  %% �����ʾ %%%%
  fprintf('Iterations:%3d  DiamMAE:%5.4e  viscIteration:%3d  ViscMAE:%5.4e.\n',Loop2_Cnt,DiamMAE,Loop1_Cnt,ViscMAE);
  if Loop1_Cnt==50
    break;
  end
end
%����Ӧ��������¼������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����Ӧ�����Ȳ����Ż����۲��� %%%%
%Ed���ܾ��仯���죻Ev�����ٱ仯����
umDiam=Diam*1e3;
%���ټ���
Vel=4.*1e3.*MeanFlow./(60.*pi.*umDiam.*umDiam);   %Ѫ������mm/s
%Ref. Bettina 2009 (��������
Ev=sqrt(mean(4*(DataArray(:,8)-Vel).^2./(DataArray(:,8)+Vel).^2));

%���ڷ�����������Ev�ͷ�
if AdapLoopOutType==2
  Ev=Ev*1.5;
else
  Ev=Ev*2;
end
y=Ev;
fprintf('Ev=%5.4e\n',y);
end


