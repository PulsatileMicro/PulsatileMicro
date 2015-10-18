function Sp = CalcSp( P )
%SP Summary of this function goes here
%   Detailed explanation goes here
taup=100-86*exp(-5000*(log(log(P)))^5.4);
Sp=-log(taup);
end

