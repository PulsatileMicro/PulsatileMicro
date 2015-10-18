function DampFactorName = GenDampFactorName(DampFactor,ViscRatio,ERatio)
global DampInit DampVisc DampE
switch DampFactor
  case DampInit
    DampFactorName=['DampInit'];
  case DampVisc
    DampFactorName=['Visc' num2str(ViscRatio)];
  case DampE
    DampFactorName=['E' num2str(ERatio)];
end