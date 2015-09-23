%% ��Ͷ�ѭ�������ж�
function [ViscMAE,LoopOutType,LoopCnt]=ViscLoopTerminator(DebugVisc,LoopCnt,T,Eju)
% MAE: max absolute error
% LoopOutType:0-������1-������2-�ﵽ����δ������3-Visc����4-���Է��̴���
LoopOutType=0;
ViscMAE=1;
if T>2  % ѭ�������㹻������������
  % ���ε������������
  Error=(DebugVisc(:,T)-DebugVisc(:,T-1))./DebugVisc(:,T-1);
  ViscMAE=max(abs(Error));
  if ViscMAE<1e-2  %�ﵽ��������
    LoopCnt=0;
    LoopOutType=1;  % ����
  elseif ~(isreal(DebugVisc(:,T))) || ~isempty(find(DebugVisc(:,T)<0, 1)) %�����жϣ�ճ�Ͷȳ��ָ���/ճ�Ͷȳ��ָ�ֵ��
    LoopCnt=0;
    LoopOutType=3;  % Visc����
  elseif Eju>0  %�����жϣ����Է��̼�����ִ���
    LoopCnt=0;
    LoopOutType=4;
  else
    LoopOutType=0;
  end
else  %ѭ����������
  if ~(isreal(DebugVisc(:,T))) || ~isempty(find(DebugVisc(:,T)<0, 1)) %�����жϣ�ճ�Ͷȳ��ָ���/ճ�Ͷȳ��ָ�ֵ��
    LoopCnt=0;
    LoopOutType=3;
  elseif Eju>0  %�����жϣ����Է��̼�����ִ��󣩣������Eju�ʹ���
    LoopCnt=0;
    LoopOutType=4;
  end
end

% ��������������޵�δ����
if LoopCnt==T && LoopCnt~=0
  LoopOutType=2;
end
