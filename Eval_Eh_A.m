%% Eval Eh, Area for 1D model
% Unit: kg m s
function [A Eh]=Eval_Eh_A(Diam, E, WallTh)
Eh=E*WallTh;
A=pi*(Diam/2).^2;
end