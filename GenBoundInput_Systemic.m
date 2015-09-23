%% ����ģ������߽����������
function [Q]=GenBoundInput_Systemic(inFileName,ModelParam)
% ȫ�ֱ���������������Macro.m�У�
dt=ModelParam(1);
Period=ModelParam(2);
NumCycles=ModelParam(3);

Qpeak=499e-6; % 485mL/s
tau=0.25;
t=dt:dt:Period;
Q1=zeros(length(t),1);
for i=1:length(t)
  if t(i)<=0.265
    Q1(i)=Qpeak*sin(pi*t(i)/tau);
  elseif t(i)>0.265 && t(i)<0.35
    Q1(i)=Q1(i-1)+Qpeak*0.0022;
  else
    Q1(i)=0;
  end
end
Q=[];
for i=1:NumCycles
  Q=[Q;Q1];
end

% ��������.bcs�ļ�
fid=fopen([inFileName '.bcs'], 'w');
fwrite(fid,Q,'double',0,'l');
fclose(fid);
