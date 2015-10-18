function [Start_Ind] = Pulse_Start(Pulse)
% ¼ì²âÂö²«²¨Æðµã(Uµã)
% Ref. Chiu YC, Arand PW, Shroff SG, Feldman T, and Carroll JD. 
% Determination of pulse wave velocities with computerized algorithms. 
% Am Heart J 121: 1460-1470, 1991.
% Pulse=Pulse(round(0.2*length(Pulse)):round(0.8*length(Pulse)));
Pulse=Pulse(200:end-200);
DerivatedPulse=diff(diff(Pulse));
[maxValue Start_Ind]=max(DerivatedPulse);
% [minValue Start_Ind]=min(Pulse);
% Start_Ind=Start_Ind+round(0.2*length(Pulse))-1;
Start_Ind=Start_Ind+200-3;
end

