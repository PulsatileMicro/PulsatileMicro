%% FahraeuLindqvist Effect
% Pries, A. R., and Secomb, T. W., 2005, "Microvascular blood viscosity in vivo and the endothelial surface layer,"
% Am. J. Physiol. Heart Circ. Physiol., 289(6), pp. H2657-H2664.
function Uvivo=FL_effect(Hd,D,Dcrit)
Doff=2.4;
% Dcrit=10.5;
D50=100;
Eamp=1.1;
Ehd=1.18;
Ewidth=0.03;
Epeak=0.6;
Wmax=2.6;

if D<=Doff
    Was=0;
else
    Was=((D-Doff)./(D+D50-2.*Doff)).*Wmax;
end

if  D<=Doff
    Wpeak=0;
elseif D>Dcrit
    Wpeak=Eamp.*exp(-Ewidth.*(D-Dcrit));
else
    Wpeak=Eamp.*(D-Doff)./(Dcrit-Doff);
end

Wph=Was+Wpeak.*Epeak;
Weff=Was+Wpeak.*(1+Hd.*Ehd);
Dph=D-2.*Wph;
Deff=D-2.*Weff;

U0=220.*exp(-1.3.*Dph)+3.2-2.44.*exp(-0.06.*(Dph.^0.645));
a=(0.8+exp(-0.075.*Dph)).*(-1+1./(1+10.^(-11).*Dph.^12))+1./(1+10.^(-11).*Dph.^12);
Uvivo=1+(U0-1).*((1-Hd).^a-1)./(0.55.^a-1);
Uvivo=Uvivo.*((D./Deff).^4);
