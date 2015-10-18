function [DiamMAE,TLoopOutType,CTimeT]=AdapLoopTerminator(DebugDiam,CTimeT,TT,Eju,AccuracyType)
%MAE: max absolute error
%LoopOutType:0-������1-������2-�ﵽ����δ������3-visc����4-���Է��̴���
TLoopOutType=0;
DiamMAE=1;
if TT<3
  if ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %�����жϣ�ճ�Ͷȳ��ָ���/ճ�Ͷȳ��ָ�ֵ��
    CTimeT=0;
    TLoopOutType=3;
  elseif Eju>0
    CTimeT=0;
    TLoopOutType=4;
  end
else
  %�ж�����ѭ��
  switch num2str(AccuracyType)
    case '0'
      Error=DebugDiam(:,TT)-DebugDiam(:,TT-1);
      [ErrorNum,ErrorIndex]=sort(abs(Error),'descend');
      %���ε�����ճ�Ͷȵ����������
      DiamMAE=ErrorNum(1);
      if DiamMAE<1e-3  %�ﵽ��������
        CTimeT=0;
        TLoopOutType=1;
      elseif ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %�����жϣ�ճ�Ͷȳ��ָ���/ճ�Ͷȳ��ָ�ֵ��
        CTimeT=0;
        TLoopOutType=3;
      elseif Eju>0  %�����жϣ����Է��̼�����ִ���
        CTimeT=0;
        TLoopOutType=4;
      else
        TLoopOutType=0;
      end
    case '1'
      Error=abs(DebugDiam(:,TT)-DebugDiam(:,TT-1))./DebugDiam(:,TT-1);
      [ErrorNum,ErrorIndex]=sort(abs(Error),'descend');
      %���ε�����ճ�Ͷȵ����������
      DiamMAE=ErrorNum(1);
      if DiamMAE<5*1e-8  %�ﵽ��������
        CTimeT=0;
        TLoopOutType=1;
      elseif ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %�����жϣ�ճ�Ͷȳ��ָ���/ճ�Ͷȳ��ָ�ֵ��
        CTimeT=0;
        TLoopOutType=3;
      elseif Eju>0  %�����жϣ����Է��̼�����ִ���
        CTimeT=0;
        TLoopOutType=4;
      else
        TLoopOutType=0;
      end
    case '2'
      Error=abs(DebugDiam(:,TT)-DebugDiam(:,TT-1));
      %���ε������ƽ���ܾ��������
      DiamMAE=mean(Error);
      if DiamMAE<1e-8  %�ﵽ��������
        CTimeT=0;
        TLoopOutType=1;
      elseif ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %�����жϣ�ճ�Ͷȳ��ָ���/ճ�Ͷȳ��ָ�ֵ��
        CTimeT=0;
        TLoopOutType=3;
      elseif Eju>0  %�����жϣ����Է��̼�����ִ���
        CTimeT=0;
        TLoopOutType=4;
      else
        TLoopOutType=0;
      end
    case '3'
      Error=abs(DebugDiam(:,TT)-DebugDiam(:,TT-1))./DebugDiam(:,TT-1);
      %���ε�����ƽ���ܾ���������
      DiamMAE=mean(Error);
      if DiamMAE<1e-3  %�ﵽ��������
        CTimeT=0;
        TLoopOutType=1;
      elseif ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %�����жϣ�ճ�Ͷȳ��ָ���/ճ�Ͷȳ��ָ�ֵ��
        CTimeT=0;
        TLoopOutType=3;
      elseif Eju>0  %�����жϣ����Է��̼�����ִ���
        CTimeT=0;
        TLoopOutType=4;
      else
        TLoopOutType=0;
      end
    otherwise
      error(' Wrang setting of "AccuracyType" !');
  end
end

if CTimeT==TT && CTimeT~=0
  TLoopOutType=2;
end

end