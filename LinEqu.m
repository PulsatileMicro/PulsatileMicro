function [MeanP,Flow,Vel,DeltaP,Visc,Hd]=LinEqu(NetTypeID,DatMatrix,Boundary,DampPara)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID
global DampInit DampVisc

%% 1.���ݽ��
% 1.1 DatMatrix
% SegName=DatMatrix(:,1);
From=DatMatrix(:,2);
To=DatMatrix(:,3);
Len=DatMatrix(:,4);
Diam=DatMatrix(:,5);
% WallTh=DatMatrix(:,6);
% SegType=DatMatrix(:,7);
Visc=DatMatrix(:,8);
% E=DatMatrix(:,9);
VesNum=length(From);

% 1.2 Boundary
BoundNode=Boundary(:,1);
BoundType=Boundary(:,2);
BoundFlow=Boundary(:,3);
BoundHd=Boundary(:,4);

% 1.3 DampPara
DampFactor=DampPara(1);
ViscRatio=DampPara(2);

%% 2.�߽�Ѫ�ܷ������� %%%%
% ���ݱ߽������ķ�������
% Ѫ��Ϊ����->��߽�->�߽�ڵ���From����
% Ѫ��Ϊ����->���߽�->�߽�ڵ���To����
[From,To]=OriNodeModify(Boundary,From,To);

Qref=0.001; % TODO: Qrefû�ã����������Ҫ�ع�
% Boundary(:,4)=0.45;
BHd=BoundaryInput(Boundary,From,To,Qref);

%% 3.��λ����
Len=Len*1e-6;   %m
Diam=Diam*1e-6;   %m
Visc=Visc*1e-3;   %Pa.s

%% 4. ��������viscosity
% ���Loop_Limit��ѭ��
if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  Loop_Limit=1;
else
  Loop_Limit=50;
end
% ѭ����������
Loop_Cnt=0;
% ��ʼHd����˳��. Porder, ����. Norder, ����
Porder=1:VesNum;Norder=VesNum:1;
% ��¼ÿ��ѭ���õ���Visc
DebugVisc=zeros(VesNum,Loop_Limit);
% ��ʼѭ��
while Loop_Cnt<Loop_Limit
  Loop_Cnt=Loop_Cnt+1;
  % ���Է������ģ��
  % ��λ�����룺SI��λ�ƣ������P:mmHg, Q:nL/min
  [MeanP,DeltaP,MeanFlow,Eju]=LinEquSolver(BoundNode,BoundType,BoundFlow,VesNum,Diam,Visc,Len,From,To);
  % ��������ģ��
  [InvIndex,FromNew,ToNew,MeanFlowNew,DeltaPNew]=AdjustFlowDir(From,To,MeanFlow,DeltaP,Eju);
  % ����Hd����˳��
  [Porder,Norder,Eju]=HdCalOrder(Boundary,From,To,FromNew,ToNew,Porder,Norder,2,Eju);
  if Eju==1
    errorFlag=1;
    break;
  end
  % ����Hd
  [Hd,FlowRatio,FQe]=HdCalc_wrf(Porder,BHd,FromNew,ToNew,Diam,MeanFlowNew,0);
  
  % ����Ӧ���������治��Ҫճ�Ͷȵ���
  if NetTypeID~=Net_546_ID && NetTypeID~=Net_122_ID
    %%ճ�Ͷȷ���
    umDiam=Diam.*1e6;   %um
    for i=1:VesNum
      Visc(i)=FL_effect(Hd(i),umDiam(i),10.5)/1e3;
      switch DampFactor
        case DampVisc
          Visc(i)=Visc(i)*ViscRatio;
      end
    end
    % ѭ������ڶ��κ�(T>1)����ʼУ��Visc
    % ����SOR����У�����������ճ�Ͷ��񵴵����
    if Loop_Cnt>1
      Visc=0.5*(Visc+DebugVisc(:,Loop_Cnt-1));
    end
  end
  
  % ��¼ÿ�ε����ķ�����
  DebugVisc(:,Loop_Cnt)=Visc;
end

% ����Ѫ������
if NetTypeID==Net_913_ID
  SegType=AutoVesType(Boundary,FromNew,ToNew,Porder);
elseif NetTypeID==Egg_818_ID
  ConstSegType=zeros(VesNum,1);
  % ����Ѫ�ܶ��޷�ͨ���Զ��㷨�ж�����
  ConstSegType([216;544;547;559;561;47;57;70;80;143;681;683])=[3;ones(11,1)];
  SegType=autoVesTypeAdv(ConstSegType,FromNew,ToNew,Porder);
elseif NetTypeID==Egg_636_ID
  ConstSegType=zeros(VesNum,1);
  ConstSegType([167;168;143;144;560;561;562;212;580;581;583;410;615;449;514;452])=ones(16,1);
  ConstSegType([213;544;506;507;583;394;384;385;386;396])=2*ones(10,1);
  SegType=autoVesTypeAdv(ConstSegType,FromNew,ToNew,Porder);
  SegType(424)=3;
  SegType(428)=2;
end

%% �������
Vel=MeanFlowNew/1e12/60./(pi*0.25*Diam.^2)*1e3;
Visc=Visc*1e3;
Flow=MeanFlowNew;
DeltaP=DeltaPNew;
% ����steady state����
save('LinEquDatFile.mat','DeltaP','Visc','Flow','Hd','MeanP','Vel');
end
