%% Calcuate Vascular Resisance
function Res=CalcVesRes(Diam,Visc,Len)
Res=8.*Visc.*Len./(pi.*(Diam/2).^4);