%% Eval gamma for 1D model
% Ref. Jordi Alastruey, Pulse wave propagation in the cardiovascular
% system. Summer School on Mathematical Modeling and Numerical Simulation of the
% Cardiovascular System,A Coru~na, July 2010
function Gamma=Eval_Gamma(Diam,phi,WallTh)
Gamma = 2/3*sqrt(pi).*phi.*WallTh./(0.25*pi*Diam.^2);
end