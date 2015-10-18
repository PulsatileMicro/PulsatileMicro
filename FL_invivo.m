function Uvivo=FahraeusLindqvist1(Hd,D)
Doff=2.4;
Dcrit=10.5;
D50=100;
Eamp=1.1;
Ehd=1.18;
% Ehd=0;
% Eamp=1.23;
Ewidth=0.03;
Epeak=0.6;
Wmax=2.6;

Was=zeros(length(D),1);
Wpeak=zeros(length(D),1);
for i=1:length(D)
    if D(i)<=Doff
        Was(i)=0;
    else
        Was(i)=((D(i)-Doff)./(D(i)+D50-2.*Doff)).*Wmax;
    end
    
    if  D(i)<=Doff
        Wpeak(i)=0;
    elseif D(i)>Dcrit
        Wpeak(i)=Eamp.*exp(-Ewidth.*(D(i)-Dcrit));
    else
        Wpeak(i)=Eamp.*(D(i)-Doff)./(Dcrit-Doff);
    end
end

Wph=Was+Wpeak.*Epeak;
Weff=Was+Wpeak.*(1+Hd.*Ehd);
Dph=D-2.*Wph;
Deff=D-2.*Weff;

U0=220.*exp(-1.3.*Dph)+3.2-2.44.*exp(-0.06.*(Dph.^0.645));
a=(0.8+exp(-0.075.*Dph)).*(-1+1./(1+10.^(-11).*Dph.^12))+1./(1+10.^(-11).*Dph.^12);
Uvivo=1+(U0-1).*((1-Hd).^a-1)./(0.55.^a-1);
Uvivo=Uvivo.*((D./Deff).^4);
