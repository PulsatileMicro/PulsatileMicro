function Sm = CalcSm( Q,Qref,Jm )
%SM Summary of this function goes here
%   Detailed explanation goes here
Sm=log(1+Jm/(Q+Qref));
end

