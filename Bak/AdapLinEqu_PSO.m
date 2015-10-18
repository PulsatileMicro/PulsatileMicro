function [AdapPara,AllEv]=AdapLinEqu(NetTypeID,DampFactor,VesNum,SegType,From,To,Diam,Len,Visc,Boundary,AdapType,FuncPara,DataArray)
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
OrgDiam=Diam;
OrgVisc=Visc;

% ����Ӧ��������ֵ
Qref=0.001;    %������������������ģ��
PO2ref=94.4;   %����ѹ���ղ���������ģ��
Tauref=0.095;  %����������������������ģ��
Jo=7142.8;     %�����źż������������ģ��
Lref=24530;    %Ѫ�ܳ���˥������������ģ��
kp=0.71043;    %ѹ���ź������Ȳ���������Ӧ�ź�
km=0.6655;     %�����ź������Ȳ���������Ӧ�ź�
kc=2.069;      %�����ź������Ȳ���������Ӧ�ź�
ks=1.6687;     %�����ź������Ȳ���������Ӧ�ź�
t=0.1;
Mo=1000;

if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  % if 0
  % ���������Ӧ���ݣ���ôֻ��һ��
  Loop1_Visc_Num=1;
  Loop2_Adap_Num=1;
  Loop3_PSO_Num=1;
  % ��ȡFuncPara���������Ϊ�߽�
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo,FuncPara);
else
  % ���ѭ������ TODO(panqing):�������ж�
  Loop1_Visc_Num=5;
  Loop2_Adap_Num=5;
  Loop3_PSO_Num=5;
  % ��FuncPara
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,Jo);
end

% PSO����
nPSO=5;
nParDim=4;
AdapPara=zeros(nParDim,nPSO,Loop3_PSO_Num);   % ��������
v=zeros(nParDim,nPSO,Loop3_PSO_Num);          % �����ٶ�
for i=1:nPSO
  AdapPara(:,i,1)=[1+rand-rand 1+rand-rand 1+rand-rand 1+rand-rand];
end
c1=4; % ��֪
c2=2; % ���
r=1;  % Լ������
w=0.729;  % ����Ȩ��
AllEv=100*ones(nPSO,Loop3_PSO_Num);   % ����Ev��Ŀ�꺯��������ֵ��Ϊ�ܴ��100

% ��¼ÿ��ѭ���õ��Ĳ���
DebugVisc=zeros(VesNum,Loop1_Visc_Num);
DebugHd=zeros(VesNum,Loop1_Visc_Num);
DebugP=zeros(VesNum,Loop1_Visc_Num);
DebugFlow=zeros(VesNum,Loop1_Visc_Num);

% ��ʼHd����˳��. Porder, ����. Norder, ����
Porder=1:VesNum;
Norder=VesNum:-1:1;

Loop3_Cnt=0;
while Loop3_Cnt<Loop3_PSO_Num
  Loop3_Cnt=Loop3_Cnt+1;
  for j=1:nPSO
    kc=AdapPara(1,j,Loop3_Cnt);
    km=AdapPara(2,j,Loop3_Cnt);
    kp=AdapPara(3,j,Loop3_Cnt);
    ks=AdapPara(4,j,Loop3_Cnt);
    Visc=OrgVisc;
    Loop2_Cnt=0;
    while Loop2_Cnt<Loop2_Adap_Num % ѭ��2������Ӧ���ؼ���
      Loop2_Cnt=Loop2_Cnt+1;
      Loop1_Cnt=0;
      Loop1_Visc_Num=50;
      
      % ���ݼ�¼����
      DebugConvergence=zeros(Loop1_Visc_Num-2,1);
      DebugVisc=zeros(VesNum,Loop1_Visc_Num);
      DebugHd=zeros(VesNum,Loop1_Visc_Num);  %��¼ѭ��2�е�Hd�仯
      DebugMeanFlow=zeros(VesNum,Loop1_Visc_Num);
      DebugDeltaP=zeros(VesNum,Loop1_Visc_Num);
      DebugMeanP=zeros(VesNum,Loop1_Visc_Num);
      DebugFQe=zeros(VesNum,Loop1_Visc_Num);
      
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
      [Stau,Sp,Tau]=SCounter(DeltaP,Diam,Len,MeanP,Tauref,Eju);
      %�ܾ�����ģ��
      Stot=Stau+kp.*Sp+km.*(Sm+kc.*Sc)-ks;
      Diam=Diam+Stot.*Diam.*t;
      
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
      %AccuracyType�����ȿ��ƣ���0-��������1-��������2-ƽ���������
      AccuracyType=2;
      [DiamMAE,TLoopOutType,Loop1_Visc_Num]=AdapLoopTerminator(DebugDiam,Loop1_Visc_Num,Loop2_Cnt,Eju,AccuracyType);
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
    Ed=sqrt(mean((DataArray(:,5)-Diam.*1e3).^2./(DataArray(:,5).^2)));
    AllEv(j,Loop3_Cnt)=Ed;
    
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
    v(:,k,Loop3_Cnt+1)=w.*v(:,k,Loop3_Cnt)+c1*rand*(Pi(:,k)-AdapPara(:,k,Loop3_Cnt))+c2*rand*(Pg-AdapPara(:,k,Loop3_Cnt));
    AdapPara(:,k,Loop3_Cnt+1)=AdapPara(:,k,Loop3_Cnt)+r*v(:,k,Loop3_Cnt+1);
  end
end

end
