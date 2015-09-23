%% 黏滞度循环结束判断
function [ViscMAE,LoopOutType,LoopCnt]=ViscLoopTerminator(DebugVisc,LoopCnt,T,Eju)
% MAE: max absolute error
% LoopOutType:0-正常；1-收敛；2-达到步长未收敛；3-Visc错误；4-线性方程错误
LoopOutType=0;
ViscMAE=1;
if T>2  % 循环次数足够，进行误差计算
  % 两次迭代间的相对误差
  Error=(DebugVisc(:,T)-DebugVisc(:,T-1))./DebugVisc(:,T-1);
  ViscMAE=max(abs(Error));
  if ViscMAE<1e-2  %达到收敛精度
    LoopCnt=0;
    LoopOutType=1;  % 收敛
  elseif ~(isreal(DebugVisc(:,T))) || ~isempty(find(DebugVisc(:,T)<0, 1)) %错误判断（粘滞度出现复数/粘滞度出现负值）
    LoopCnt=0;
    LoopOutType=3;  % Visc出错
  elseif Eju>0  %错误判断（线性方程计算出现错误）
    LoopCnt=0;
    LoopOutType=4;
  else
    LoopOutType=0;
  end
else  %循环次数不够
  if ~(isreal(DebugVisc(:,T))) || ~isempty(find(DebugVisc(:,T)<0, 1)) %错误判断（粘滞度出现复数/粘滞度出现负值）
    LoopCnt=0;
    LoopOutType=3;
  elseif Eju>0  %错误判断（线性方程计算出现错误），传入的Eju就错误
    LoopCnt=0;
    LoopOutType=4;
  end
end

% 到达迭代次数上限但未收敛
if LoopCnt==T && LoopCnt~=0
  LoopOutType=2;
end
