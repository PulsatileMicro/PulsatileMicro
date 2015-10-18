function [DiamMAE,TLoopOutType,CTimeT]=AdapLoopTerminator(DebugDiam,CTimeT,TT,Eju,AccuracyType)
%MAE: max absolute error
%LoopOutType:0-正常；1-收敛；2-达到步长未收敛；3-visc错误；4-线性方程错误
TLoopOutType=0;
DiamMAE=1;
if TT<3
  if ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %错误判断（粘滞度出现复数/粘滞度出现负值）
    CTimeT=0;
    TLoopOutType=3;
  elseif Eju>0
    CTimeT=0;
    TLoopOutType=4;
  end
else
  %判断跳出循环
  switch num2str(AccuracyType)
    case '0'
      Error=DebugDiam(:,TT)-DebugDiam(:,TT-1);
      [ErrorNum,ErrorIndex]=sort(abs(Error),'descend');
      %两次迭代间粘滞度的最大绝对误差
      DiamMAE=ErrorNum(1);
      if DiamMAE<1e-3  %达到收敛精度
        CTimeT=0;
        TLoopOutType=1;
      elseif ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %错误判断（粘滞度出现复数/粘滞度出现负值）
        CTimeT=0;
        TLoopOutType=3;
      elseif Eju>0  %错误判断（线性方程计算出现错误）
        CTimeT=0;
        TLoopOutType=4;
      else
        TLoopOutType=0;
      end
    case '1'
      Error=abs(DebugDiam(:,TT)-DebugDiam(:,TT-1))./DebugDiam(:,TT-1);
      [ErrorNum,ErrorIndex]=sort(abs(Error),'descend');
      %两次迭代间粘滞度的最大相对误差
      DiamMAE=ErrorNum(1);
      if DiamMAE<5*1e-8  %达到收敛精度
        CTimeT=0;
        TLoopOutType=1;
      elseif ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %错误判断（粘滞度出现复数/粘滞度出现负值）
        CTimeT=0;
        TLoopOutType=3;
      elseif Eju>0  %错误判断（线性方程计算出现错误）
        CTimeT=0;
        TLoopOutType=4;
      else
        TLoopOutType=0;
      end
    case '2'
      Error=abs(DebugDiam(:,TT)-DebugDiam(:,TT-1));
      %两次迭代间的平均管径绝对误差
      DiamMAE=mean(Error);
      if DiamMAE<1e-8  %达到收敛精度
        CTimeT=0;
        TLoopOutType=1;
      elseif ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %错误判断（粘滞度出现复数/粘滞度出现负值）
        CTimeT=0;
        TLoopOutType=3;
      elseif Eju>0  %错误判断（线性方程计算出现错误）
        CTimeT=0;
        TLoopOutType=4;
      else
        TLoopOutType=0;
      end
    case '3'
      Error=abs(DebugDiam(:,TT)-DebugDiam(:,TT-1))./DebugDiam(:,TT-1);
      %两次迭代间平均管径的相对误差
      DiamMAE=mean(Error);
      if DiamMAE<1e-3  %达到收敛精度
        CTimeT=0;
        TLoopOutType=1;
      elseif ~(isreal(DebugDiam(:,TT))) || ~isempty(find(DebugDiam(:,TT)<0, 1)) %错误判断（粘滞度出现复数/粘滞度出现负值）
        CTimeT=0;
        TLoopOutType=3;
      elseif Eju>0  %错误判断（线性方程计算出现错误）
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