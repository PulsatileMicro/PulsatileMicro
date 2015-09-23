function SolverName = GetSolverName(Solver)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP Sparse_SS
switch Solver
  case ONED_EXP
    SolverName='1DEXP';
  case ONED_IMP
    SolverName='1DIMP';
  case RLC_EXP
    SolverName='RLCEXP';
  case RLC_IMP
    SolverName='RLCIMP';
  case SS
    SolverName='SS';
  case WOM_1
    SolverName='WOM_1';
  case WOM_2
    SolverName='WOM_2';
  case RC_EXP
    SolverName='RCEXP';
  case RC_IMP
    SolverName='RCIMP';
  case Sparse_SS
    SolverName='SparSS';
end

