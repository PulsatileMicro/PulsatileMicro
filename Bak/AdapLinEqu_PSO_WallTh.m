function [AdapPara,AllEv,AllWallTh,AllDiam,Visc,Sm,Sc]=AdapLinEqu(NetTypeID,DampFactor,VesNum,SegType,From,To,Diam,Len,Visc,WallTh,Boundary,AdapType,FuncPara,DataArray)
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
WallTh=WallTh*1e-3;
% ����ԭ����Diam��������Edʹ��
OrgDiam=Diam;
OrgVisc=Visc;
OrgWallTh=WallTh;

% ����Ӧ��������ֵ
Qref=0.001;    %������������������ģ��
PO2ref=94.4;   %����ѹ���ղ���������ģ��
tauref=0.5598;  %����������������������ģ��
Jo=6618;     %�����źż������������ģ��
Lref=14292;    %Ѫ�ܳ���˥������������ģ��
t=0.01;
Mo=1000;

if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  % if 0
  % ���������Ӧ���ݣ���ôֻ��һ��
  Loop1_ViscNum_Init=1;
  Loop2_Adap_Num=1;
  Loop3_PSO_Num=1;
  % ��ȡFuncPara���������Ϊ�߽�
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo,FuncPara);
else
  % ���ѭ������ TODO(panqing):�������ж�
  Loop1_ViscNum_Init=50;
  Loop2_Adap_Num_Init=2000;
  Loop3_PSO_Num=1;
  % ��FuncPara
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo);
end

% PSO����
nPSO=1;
nParDim=7;
AdapPara=zeros(nParDim,nPSO,Loop3_PSO_Num);   % ��������
v=zeros(nParDim,nPSO,Loop3_PSO_Num);          % �����ٶ�
AllDiam=zeros(nPSO,VesNum);
AllWallTh=zeros(nPSO,VesNum);
InitAdapPara=[1.66,0.955,-0.374,3.077,0.0177,0.114,0.609];
for i=1:nPSO
%   AdapPara(:,i,1)=InitAdapPara+0.05*rand(1,nParDim)-0.05*rand(1,nParDim);
  AdapPara(:,i,1)=InitAdapPara;
end
c1=2; % ��֪
c2=2; % ���
r=1;  % Լ������
ww=0.729;  % ����Ȩ��
AllEv=100*ones(nPSO,Loop3_PSO_Num);   % ����Ev��Ŀ�꺯��������ֵ��Ϊ�ܴ��100

% ��¼ÿ��ѭ���õ��Ĳ���
DebugVisc=zeros(VesNum,Loop1_ViscNum_Init);
DebugHd=zeros(VesNum,Loop1_ViscNum_Init);
DebugP=zeros(VesNum,Loop1_ViscNum_Init);
DebugFlow=zeros(VesNum,Loop1_ViscNum_Init);

Loop3_Cnt=0;
% PSO����
while Loop3_Cnt<Loop3_PSO_Num
  Loop3_Cnt=Loop3_Cnt+1;
  Diam=OrgDiam;
  Visc=OrgVisc;
  WallTh=OrgWallTh;
  Loop2_Cnt=0;
  for j=1:nPSO
    kc=AdapPara(1,j,Loop3_Cnt);
    kmd=AdapPara(2,j,Loop3_Cnt);
    kmg=AdapPara(3,j,Loop3_Cnt);
    ksd=AdapPara(4,j,Loop3_Cnt);
    ksg=AdapPara(5,j,Loop3_Cnt);
    kwt=AdapPara(6,j,Loop3_Cnt);
    kwo=AdapPara(7,j,Loop3_Cnt);
      
    % ÿ������Ӧ����ǰ���轫�����ָ�ΪĬ��ֵ
    Loop2_Cnt=0;
    Loop2_Adap_Num=Loop2_Adap_Num_Init;
    Diam=OrgDiam;
    Visc=OrgVisc;
    WallTh=OrgWallTh;
    while Loop2_Cnt<Loop2_Adap_Num % ѭ��2������Ӧ���ؼ���
      Loop2_Cnt=Loop2_Cnt+1;
      Loop1_Cnt=0;
      Loop1_Visc_Num=Loop1_ViscNum_Init;
      
      % ���ݼ�¼����
      DebugConvergence=zeros(Loop1_Visc_Num-2,1);
      DebugVisc=zeros(VesNum,Loop1_Visc_Num);
      DebugHd=zeros(VesNum,Loop1_Visc_Num);  %��¼ѭ��2�е�Hd�仯
      DebugMeanFlow=zeros(VesNum,Loop1_Visc_Num);
      DebugDeltaP=zeros(VesNum,Loop1_Visc_Num);
      DebugMeanP=zeros(VesNum,Loop1_Visc_Num);
      DebugFQe=zeros(VesNum,Loop1_Visc_Num);
      
      LoopOutFlag=0;
      % ��ʼHd����˳��. Porder, ����. Norder, ����
      Porder=1:VesNum;
      Norder=VesNum:-1:1;
%       Visc=OrgVisc;
      while Loop1_Cnt<Loop1_Visc_Num   % ѭ��1��Visc����
        Loop1_Cnt=Loop1_Cnt+1;
        % ���Է������ģ��
        [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam/1e3,Visc,Len/1e3,From,To);
        % ��������ģ��
        [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju);
        % ����Hd����˳��
        [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,1,Eju);
        if Eju==1
          break;
        end
        % ����Hd
        [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam/1e3,MeanFlowNew,Eju);
        
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
      
      %PO2����ģ��
      [PO2,SO2in,SO2mid,SO2out]=PO2Counter(Porder,BSO2in,FromNew,ToNew,Hd,Len,MeanFlowNew,Eju);
      %Sm����ģ��
      [Sm,Jm]=SmCounter(Porder,BJm,FromNew,ToNew,Len,MeanFlowNew,PO2,Qref,PO2ref,t,Mo,Eju);
      %Sc����ģ��
      [Sc,Jc]=ScCounter(Norder,BJc,FromNew,ToNew,Len,Sm,Jo,Lref,Eju);
      %Stau��Sp����ģ��
      % TODO: ��Ҫ�޸�ΪHypertension�����еķ���
      [Stau,Sp,Tau]=SCounter(DeltaP,Diam,Len,MeanP,tauref,Eju);
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
      Oref=32050;  %dyn/cm2
      wref=0.804;  %um
      ee=0.01;
      % �����㣺Rtd^2+Rod^2=1, Rtg^2+Rog^2=1
      Rtd=0.8;Rod=0;Rtg=0;Rog=0.6;Rw=0.2;
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
      DebugStau(:,Loop2_Cnt)=Stau;
      DebugSp(:,Loop2_Cnt)=Sp;
      DebugSm(:,Loop2_Cnt)=Sm;
      DebugSc(:,Loop2_Cnt)=Sc;
      DebugDiam(:,Loop2_Cnt)=Diam*1e3;
      DebugFlow(:,Loop2_Cnt)=MeanFlowNew;
      DebugP(:,Loop2_Cnt)=MeanP;
      DebugPO2(:,Loop2_Cnt)=PO2;
      DebugTau(:,Loop2_Cnt)=Tau;
      DebugSO2in(:,Loop2_Cnt)=SO2in;
      DebugSO2out(:,Loop2_Cnt)=SO2out;
      DebugTHd(:,Loop2_Cnt)=Hd;
      
      %% ѭ������ģ�� %%%%
      %MAE: max absolute error
      %TLoopOutType:0-������1-������2-�ﵽ����δ������3-visc����4-���Է��̴���
      %AccuracyType�����ȿ��ƣ���0-��������1-��������2-ƽ��������3-ƽ���ܾ�������
      AccuracyType=3;
      [DiamMAE,TLoopOutType,Loop2_Adap_Num]=AdapLoopTerminator(DebugDiam,Loop2_Adap_Num,Loop2_Cnt,Eju,AccuracyType);
      %�쳣���������������
      if TLoopOutType>2
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
      fprintf('PSOIter=%d, PSOParticle=%d, AdapIter=%d, DiamMAE:%5.4e,ViscIter=%d,ViscMAE:%5.4e.\n',Loop3_Cnt,j,Loop2_Cnt,DiamMAE,Loop1_Cnt,ViscMAE);
      if Loop1OutType~=1
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
    if isempty(DataArray)
      Ed=sqrt(mean((OrgDiam-Diam).^2./(OrgDiam.^2)));
    else
      Ev=sqrt(mean(4*(DataArray(:,8)-Vel).^2./(DataArray(:,8)+Vel).^2));
      Ed=sqrt(mean(4*(DataArray(:,5)-Diam*1e3).^2./(DataArray(:,5).^2)));
    end
    if LoopOutFlag==1
      %       Ev=Ev*10;
      Ed=Ed*10;
    end
    AllEv(j,Loop3_Cnt)=Ed;
    AllDiam(j,:)=Diam*1e3;
    AllWallTh(j,:)=WallTh*1e3;
    
    %     %% ����Ӧ������� %%%%
    %     AdaptationResult=[SegName DataArray(:,5) umDiam MeanP Tau FlowMean Hd PO2 Stot Stau Sp Sm Sc DeltaP];
    %     Parameters=[TLoopOutType AccuracyType DiamMAE Qref PO2ref Tauref Jo Lref kp km kc ks Ev Ed];
    %     AdaptData=[Parameters;AdaptationResult];
    %     save('AdaptData.mat','AdaptData');
    %     clear Vel umDiam;
  end
  % ������Ӧ�Ⱥ�����������һ�ε���������
  % 1. Ѱ����ʷ�����е�ȫ������ֵ
  [Value1,Ind1]=min(AllEv);
  [Value2,Ind2]=min(Value1);
  Pg=AdapPara(:,Ind1(Ind2),Ind2);
  % 2. Ѱ����ʷ�����У�ÿ�����ӵ�����ֵ
  for k=1:nPSO
    [Value,Ind]=min(AllEv(k,:));   % IndΪ��k�����ӵ���ʷ����ֵ�����
    Pi(:,k)=AdapPara(:,k,Ind);  % Pi(:,k)Ϊ��k�����ӵ���ʷ����ֵ
    v(:,k,Loop3_Cnt+1)=ww.*v(:,k,Loop3_Cnt)+c1*rand*(Pi(:,k)-AdapPara(:,k,Loop3_Cnt))+c2*rand*(Pg-AdapPara(:,k,Loop3_Cnt));
    AdapPara(:,k,Loop3_Cnt+1)=AdapPara(:,k,Loop3_Cnt)+r*v(:,k,Loop3_Cnt+1);
  end
end

