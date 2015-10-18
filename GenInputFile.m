%% 生成输入文件
function fileName=GenInputFile(VesType,VesParam,BCType,Bifur,BCVal,ModelParam,DampFactorName,dt,nCycle,Period,len_ratio,mass_ratio)
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP Sparse_SS
global NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE IN_BOUND_TYPE RANDOM_PHASE VEL_PROFILE USE_OPTBOUND BIFURLOSS
Len=VesParam(1,:);
Diam=VesParam(2,:);
WallTh=VesParam(3,:);
E=VesParam(4,:);
Phi=VesParam(5,:);
alpha=VesParam(6,:);
visc=VesParam(7,:);
Hd=VesParam(8,:);
SegName=VesParam(9,:);
From=VesParam(10,:);
To=VesParam(11,:);
Vel=VesParam(12,:);
Gamma=VesParam(13,:);
GammaII=VesParam(14,:);
q=VesParam(15,1);
L=VesParam(16,1);
scale_lamda=VesParam(17,1);
scale_u0=VesParam(18,1);
scale_r0=VesParam(19,1);
NumHisPt=VesParam(20,1);
NoDim=ModelParam(NONDIM);
UpdVisc=ModelParam(VISC_UPD);
Riem=ModelParam(RIEM);
Binout=ModelParam(BINOUT);
ViscWall=ModelParam(VISCOELAS);
CFL_Flag=ModelParam(CFL);
StepLapse=ModelParam(STEP_LAPSE);
ODESolver=ModelParam(ODESOLVER);
InputWave=ModelParam(INPUT_WAVE);
RelTol=1e-6;
AbsTol=1e-6;
BifurLoss=ModelParam(BIFURLOSS);

VesNum=length(Len);
fileName=[VesType '.in'];
fid=fopen(fileName, 'w');

if NoDim
  fprintf(fid, '19 parameter list\n');
else
  fprintf(fid, '16 parameter list\n');
end
fprintf(fid, '%E\tDT\n', dt);
if ModelParam(ODESOLVER)==SS || ModelParam(ODESOLVER)==Sparse_SS % Steady state model
  fprintf(fid, '%d\tNSTEPS\n', 1);
else
  fprintf(fid, '%d\tNSTEPS\n', uint32(Period/dt*nCycle));
end
if ModelParam(ODESOLVER)~=0
  fprintf(fid, '%lu\tHISSTEP\n', 1);
else
  fprintf(fid, '%lu\tHISSTEP\n', uint32(1/dt/1000));
end
if ModelParam(ODESOLVER)==RLC_IMP || ModelParam(ODESOLVER)==RLC_EXP ||...
    ModelParam(ODESOLVER)==RC_IMP || ModelParam(ODESOLVER)==RC_EXP
  fprintf(fid, '%E\tRho (kg/m3)\n',1.050*mass_ratio*len_ratio^(-3));
else
  fprintf(fid, '%E\tRho (kg/m3)\n',1050);
end
fprintf(fid, '%f\tAlpha (velocity profile param in the momentum eq) (for Alpha=1 no wall viscous effects considered)\n', alpha(1));
fprintf(fid, '%f\tpinf(Pa)\n', 0);
fprintf(fid, '%d\tUpdVisc\n', UpdVisc);
fprintf(fid, '%d\tRiemann\n', Riem);
fprintf(fid, '%d\tBinout\n', Binout);
fprintf(fid, '%E\tRELTOL\n',RelTol);
fprintf(fid, '%E\tABSTOL\n',AbsTol);
fprintf(fid, '%d\tSHOWCFL\n',CFL_Flag);
fprintf(fid, '%d\tSTEPLAPSE\n',StepLapse);
fprintf(fid, '%d\tODESolver\n',ODESolver);
fprintf(fid, '%d\tBifurLoss\n',BifurLoss);
if Gamma(1)~=0
  fprintf(fid, '1\tGAMMA\n');
elseif GammaII(1)~=0
  fprintf(fid, '2\tGAMMA\n');
else
  fprintf(fid, '0\tGAMMA\n');
end
if NoDim
  fprintf(fid, '%E\tSCALE_LAMDA\n',scale_lamda);
  fprintf(fid, '%E\tSCALE_U0\n',scale_u0);
  fprintf(fid, '%E\tSCALE_R0\n',scale_r0);
end

%%%% Mesh section %%%%
fprintf(fid, 'Mesh  -- expansion order -- quadrature order Ndomains = %d\n', VesNum);
for i=1:VesNum
  [A(i) Eh(i)]=Eval_Eh_A(Diam(i), E(i), WallTh(i));
  if NoDim
    % Non-dimensionalization
    A(i)=A(i)/scale_r0/scale_r0;
    Len(i)=Len(i)/scale_lamda;
    Eh(i)=Eh(i)/1050/scale_u0/scale_u0/scale_r0;
  end
  if Gamma(1) ~= 0
    fprintf(fid, '1\tnel Eh Area Visc Hd Gamma\n');
  elseif GammaII(1) ~= 0
    fprintf(fid, '1\tnel Eh Area Visc Hd GamII\n');
  else
    fprintf(fid, '1\tnel Eh Area Visc Hd\n');
  end
  fprintf(fid, '0.0\t %.12f \t%d\t%d\t# x_lower x_upper L q\n', Len(i), L, q);
  fprintf(fid, 'Ao = %.16E\n', A(i));
  fprintf(fid, 'Visc = %.12f\n', visc(i));
  fprintf(fid, 'Hd = %.12f\n', Hd(i));
  if Gamma(1) ~= 0
    if NoDim
      fprintf(fid, 'Gamma = %f\n', Gamma(i)/1050/scale_u0/scale_r0/scale_lamda);
    else
      fprintf(fid, 'Gamma = %f\n', Gamma(i));
    end
  elseif GammaII(1) ~= 0
    if NoDim
      fprintf(fid, 'GamII = %f\n', GammaII(i)/1050/scale_u0/scale_r0/scale_lamda);
    else
      fprintf(fid, 'GamII = %f\n', GammaII(i));
    end
  end
  fprintf(fid, 'Eh = %.12f\n', Eh(i));
end
fprintf(fid, 'Boundary conditons\n');
for i=1:VesNum
  switch BCType(i,1)
    case 'u'
      fprintf(fid, 'u 0     # A lhs boundary Domain %d\n', i);
      fprintf(fid, '   a = %.20E\n', A(i));
%       fprintf(fid, '   a = Ao\n');
      fprintf(fid, 'u %d     # U lhs boundary\n', Bifur(i,2));
      fprintf(fid, '   u = %.12E\n', -1);   % Velocity in mm/s
      fprintf(fid, 'h 0     # Hd lhs boundary\n');
      fprintf(fid, '   h = %f\n', Bifur(i,1));
    case 'q'
      fprintf(fid, 'q 0     # A lhs boundary Domain %d\n',i);
      fprintf(fid, '   a = %.20E\n', A(i));
      if ModelParam(ODESOLVER)~=SS && ModelParam(ODESOLVER)~=Sparse_SS
        fprintf(fid, 'q %d     # Q lhs boundary\n', Bifur(i,2));
      else
        fprintf(fid, 'q %.20E     # Q lhs boundary\n', BCVal(i));
      end
      fprintf(fid, '   q = %.20E\n', -1);  % Flow in mm^3/s
      fprintf(fid, 'h 0     # Hd lhs boundary\n');
      fprintf(fid, '   h = %f\n', Bifur(i,1));
    case 'p'
      fprintf(fid, 'p %d     # Q lhs boundary\n', Bifur(i,2));
      fprintf(fid, '   q = %.20E\n', -1);
      fprintf(fid, 'p %d     # A lhs boundary Domain %d\n',Bifur(i,2), i);
      fprintf(fid, '   u = %.20E\n', Vel(i)/1e3);
%       fprintf(fid, '   u = %.20E\n', 0);
      fprintf(fid, 'h 0     # Hd lhs boundary\n');
      fprintf(fid, '   h = %f\n', Bifur(i,1));
    case 'B'
      fprintf(fid, 'B  %d  %d  # A rhs boundary Bifur with %d %d  Domain %d\n', Bifur(i,1), Bifur(i,2), Bifur(i,1), Bifur(i,2), i);
      fprintf(fid, 'B  %d  %d  # U rhs boundary Bifur with %d %d\n', Bifur(i,1), Bifur(i,2), Bifur(i,1), Bifur(i,2));
      fprintf(fid, 'B  %d  %d  # Hd rhs boundary Bifur with %d %d\n', Bifur(i,1), Bifur(i,2), Bifur(i,1), Bifur(i,2));
    case 'C'
      fprintf(fid, 'C  %d  %d  # A rhs boundary Bifur with %d %d Domain %d\n', Bifur(i,1), Bifur(i,2), Bifur(i,1), Bifur(i,2), i);
      fprintf(fid, 'C  %d  %d  # U rhs boundary Bifur with %d %d\n', Bifur(i,1), Bifur(i,2), Bifur(i,1), Bifur(i,2));
      fprintf(fid, 'C  %d  %d  # Hd rhs boundary Bifur with %d %d\n', Bifur(i,1), Bifur(i,2), Bifur(i,1), Bifur(i,2));
    case 'J'
      fprintf(fid, 'J  %d  %d  # A rhs boundary Bifur with %d %d Domain %d\n', Bifur(i,1), Bifur(i,2), Bifur(i,1), Bifur(i,2), i);
      fprintf(fid, 'J  %d  %d  # U rhs boundary Bifur with %d %d\n', Bifur(i,1), Bifur(i,2), Bifur(i,1), Bifur(i,2));
      fprintf(fid, 'J  %d  %d  # Hd rhs boundary Bifur with %d %d\n', Bifur(i,1), Bifur(i,2), Bifur(i,1), Bifur(i,2));
    case 'T'
      fprintf(fid, 'T  %f\n', BCVal(i));
      fprintf(fid, 'T  %f\n', BCVal(i));
      fprintf(fid, 'H  0\n');
    case 't'
      fprintf(fid, 't  %f\n', BCVal(i));
      fprintf(fid, 't  %f\n', BCVal(i));
      fprintf(fid, 'H  0\n');
  end
  switch BCType(i,2)
    case 'B'
      fprintf(fid, 'B  %d  %d  # A rhs boundary Bifur with %d %d\n', Bifur(i,3), Bifur(i,4), Bifur(i,3), Bifur(i,4));
      fprintf(fid, 'B  %d  %d  # U rhs boundary Bifur with %d %d\n', Bifur(i,3), Bifur(i,4), Bifur(i,3), Bifur(i,4));
      fprintf(fid, 'B  %d  %d  # Hd rhs boundary Bifur with %d %d\n', Bifur(i,3), Bifur(i,4), Bifur(i,3), Bifur(i,4));
    case 'C'
      fprintf(fid, 'C  %d  %d  # A rhs boundary Bifur with %d %d\n', Bifur(i,3), Bifur(i,4), Bifur(i,3), Bifur(i,4));
      fprintf(fid, 'C  %d  %d  # U rhs boundary Bifur with %d %d\n', Bifur(i,3), Bifur(i,4), Bifur(i,3), Bifur(i,4));
      fprintf(fid, 'C  %d  %d  # Hd rhs boundary Bifur with %d %d\n', Bifur(i,3), Bifur(i,4), Bifur(i,3), Bifur(i,4));
    case 'J'
      fprintf(fid, 'J  %d  %d  # A rhs boundary Bifur with %d %d\n', Bifur(i,3), Bifur(i,4), Bifur(i,3), Bifur(i,4));
      fprintf(fid, 'J  %d  %d  # U rhs boundary Bifur with %d %d\n', Bifur(i,3), Bifur(i,4), Bifur(i,3), Bifur(i,4));
      fprintf(fid, 'J  %d  %d  # Hd rhs boundary Bifur with %d %d\n', Bifur(i,3), Bifur(i,4), Bifur(i,3), Bifur(i,4));
    case 'R'
      if NoDim
        BCVal(i)=BCVal(i)/1050/scale_u0*scale_r0*scale_r0;
      end
      fprintf(fid, 'R  %.16E\n', BCVal(i));
      fprintf(fid, 'R  %.16E\n', BCVal(i));
      fprintf(fid, 'H  0\n');
    case 'W'
      fprintf(fid, 'W %.16E\n', Bifur(i,3));
      fprintf(fid, 'W %.16E\n', Bifur(i,4));
      fprintf(fid, 'H  0\n');
    case 'T'
      fprintf(fid, 'T  %f\n', BCVal(i));
      fprintf(fid, 'T  %f\n', BCVal(i));
      fprintf(fid, 'H  0\n');
    case 't'
      fprintf(fid, 't  %f\n', BCVal(i));
      fprintf(fid, 't  %f\n', BCVal(i));
      fprintf(fid, 'H  0\n');
    case 'q'
      fprintf(fid, 'q 0     # A rhs boundary\n');
      fprintf(fid, '   a = %E\n', A(i));
      if ModelParam(ODESOLVER)~=SS && ModelParam(ODESOLVER)~=Sparse_SS
        fprintf(fid, 'q %d     # Q lhs boundary\n', Bifur(i,2));
      else
        fprintf(fid, 'q %.20E     # Q lhs boundary\n', BCVal(i));
      end
      fprintf(fid, '   q = %E\n', -1);
      fprintf(fid, 'h 0     # Hd rhs boundary\n');
      fprintf(fid, '   h = 0.45\n');
    case 'u'
      fprintf(fid, 'u 0      # A rhs boundary\n');
      fprintf(fid, '   a = %.20E\n', A(i));
      fprintf(fid, 'u %d     # U rhs boundary\n', Bifur(i,4));
      fprintf(fid, '   u = %.12f\n', BCVal(i));
      fprintf(fid, 'h 0     # Hd lhs boundary\n');
      fprintf(fid, '   h = %f\n', Bifur(i,1));
    case 'p'
      fprintf(fid, 'p 0      # A rhs boundary\n');
      fprintf(fid, '   a = %.20E\n', A(i));
      fprintf(fid, 'p %d     # U rhs boundary\n', Bifur(i,4));
      fprintf(fid, '   u = %.12f\n', BCVal(i));
      fprintf(fid, 'h 0     # Hd lhs boundary\n');
      fprintf(fid, '   h = %f\n', Bifur(i,1));
  end
end
fprintf(fid, 'Initial condition\n');
for i=1:VesNum
  fprintf(fid, 'a = Ao\n');
  fprintf(fid, 'u = %.8E\n', Vel(i)/1000);
%   fprintf(fid, 'u = %.6f\n', 0);
  fprintf(fid, 'h = %.6f\n', Hd(i));
end
fprintf(fid, 'History Pts\n');
fprintf(fid, '%d #Number of Domains with history points\n', VesNum);
for i=1:VesNum
  fprintf(fid, '%d %d  #Npts Domain id x[1], x[2], x[3]\n', NumHisPt, i);
  for j=1:NumHisPt
    fprintf(fid, '%.12f ', Len(i)*(j-1)/(NumHisPt-1));
  end
  fprintf(fid, '\n');
end
fclose(fid);