function Ehr_ratio = Eh_Olufsen(radius)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% radius - cm
% Eh - g/s2
k1=2e7;   % g/s2/cm1
k2=-22.53;% cm-1
k3=8.65e5;% g/s2/cm1
Ehr_ratio=k1*exp(k2*radius)+k3;
end

