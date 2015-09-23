%% Eval beta & Area for 1D model
% Unit: kg m s
function [Area beta]=Eval_beta_A(Diam, E, WallTh)
Area=pi*(Diam/2).^2;
beta=4/3*sqrt(pi)*WallTh.*E./Area;