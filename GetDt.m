%% �������������ȷ�����沽��
function dt=GetDt(ODESolverType)
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP Sparse_SS
if ODESolverType~=ONED_EXP
  % 1D��ʽ����0D
  dt=1e-3;
else
  % 1D��ʽ
  dt=1e-4;
end
