%% Wave separation analysis
function [Pp,Pn,Up,Un]=WIA(P,U,A,beta,E,h,Visc)
dataLen=length(P)-1;
rho=1050;   % kg/m3
c=sqrt(beta./2./rho).*A.^0.25;
% D=sqrt(4*A/pi);
% c=D/4.*sqrt(1.25*2*pi.*E.*h./D./Visc/0.75);
% c=35.4*D.*sqrt(1.25*2*pi.*E.*h./D./Visc);
dP=diff(P);
dU=diff(U);
c(1)=[];

dPp=(dP+rho.*c.*dU)/2;      % diff of positive pressure
dPn=(dP-rho.*c.*dU)/2;      % diff of negative pressure
dUp=(dU+dP./rho./c)/2;      % diff of positive velocity
dUn=(dU-dP./rho./c)/2;      % diff of negative velocity

Pp=zeros(dataLen,1);
Pn=zeros(dataLen,1);
Up=zeros(dataLen,1);
Un=zeros(dataLen,1);
Pref=P(1);
Uref=U(1);
for i=1:dataLen
    Pp(i)=sum(dPp(1:i))+Pref;
    Pn(i)=sum(dPn(1:i))+Pref;
    Up(i)=sum(dUp(1:i))+Uref;
    Un(i)=sum(dUn(1:i))+Uref;
end
P(1)=[];
U(1)=[];