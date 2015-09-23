%% Eval viscous modulus Phi for 1D model
% Ref. Avolio A. P., Multi-Branched model of the human arterial system [J], 
% Medical & Biological Engineering & Computing, 1980, 18(6): 709-718.
function Phi=Eval_Phi(E,Period)
Angle0=15;  % 15¶È
f=1/Period;
w=2*pi*f;
Angle=Angle0*(1-exp(-2*w));
Phi=E./w*tan(Angle/180*pi);
end