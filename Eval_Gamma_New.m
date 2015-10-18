%% Eval gamma for 1D model
% Ref. Alastruey J., et al., Pulse wave propagation in a model human arterial network: 
% Assessment of 1-D visco-elastic simulations against in vitro measurements, Journal of Biomechanics, 2011, 44(12): 2250-2258.
function Gamma=Eval_Gamma_New(Diam,phi)
Gamma = phi/2/sqrt(pi)/sqrt((0.25*pi*Diam.^2));
end