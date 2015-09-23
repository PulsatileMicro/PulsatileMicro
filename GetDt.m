%% 根据求解器类型确定仿真步长
function dt=GetDt(ODESolverType)
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP Sparse_SS
if ODESolverType~=ONED_EXP
  % 1D隐式或者0D
  dt=1e-3;
else
  % 1D显式
  dt=1e-4;
end
