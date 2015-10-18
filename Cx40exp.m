function Cex = Cx40exp(tau,Cbasal,Cpeak,tau_peak,tau_width)
% Compute the expression of Cx40
Psai=1./(1+exp(-(tau-tau_peak)./(0.14*tau_width)));
Cex=Cbasal+4*Cpeak.*Psai./((1+Psai).^2);

end

