%% Generate input file for Symetric Network
function fileName=GenSymNet(VesType, VesParam, BCType, Bifur, dt, Nstep)
Len=VesParam(1,:);
Diam=VesParam(2,:);
WallTh=VesParam(3,:);
E=VesParam(4,:);
alpha=VesParam(5,:);
visc=VesParam(6,:);
Hd=VesParam(7,:);
pinf=VesParam(8,:);
PeakFlow=VesParam(9,:);
SegName=VesParam(10,:);
From=VesParam(11,:);
To=VesParam(12,:);
Vel=VesParam(13,:);
Gamma=VesParam(14,:);
GammaII=VesParam(15,:);

VesNum=length(SegName);
fileName=[VesType '.in'];
fid=fopen(fileName, 'w');

fprintf(fid, '12 parameter list\n');
fprintf(fid, '%.12E\tDT\n', dt);
fprintf(fid, '%lu\tNSTEPS\n', uint32(Nstep));
fprintf(fid, '%lu\tIOSTEP\n', uint32(Nstep));
fprintf(fid, '%lu\tHISSTEP\n', uint32(1/dt/1000));
fprintf(fid, '1050\tRho (kg/m3)\n');
fprintf(fid, '%.12f\tAlpha (velocity profile param in the momentum eq) (for Alpha=1 no wall viscous effects considered)\n', alpha(1));
fprintf(fid, '%.12f\tpinf(Pa)\n', pinf(1));
fprintf(fid, '%.12E\tRELTOL\n',1e-6);
fprintf(fid, '%.12E\tABSTOL\n',1e-6);
fprintf(fid, '%d\tRiemann\n', 1);
fprintf(fid, '%d\tSHOWCFL\n',0);
fprintf(fid, '%d\tSTEPLAPSE\n',0);
fprintf(fid, 'Mesh  -- expansion order -- quadrature order Ndomains = %d\n', VesNum);
for i=1:VesNum
  [A(i) Eh(i)]=Eval_Eh_A(Diam(i), E(i), WallTh(i));
  if Gamma(1) ~= 0
    fprintf(fid, '1\tnel Eh Area Visc Hd Gamma\n');
  elseif GammaII(1) ~= 0
    fprintf(fid, '1\tnel Eh Area Visc Hd GamII\n');
  else
    fprintf(fid, '1\tnel Eh Area Visc Hd\n');
  end
  fprintf(fid, '0.0\t %.12f \t3\t3\t# x_lower x_upper L q\n', Len(i));
  fprintf(fid, 'Ao = %.12E\n', A(i));
  fprintf(fid, 'Visc = %.12f\n', visc(i));
  fprintf(fid, 'Hd = %.12f\n', Hd(i));
  if Gamma(1) ~= 0
    fprintf(fid, 'Gamma = %.12f\n', Gamma(i));
  elseif GammaII(1) ~= 0
    fprintf(fid, 'GamII = %.12f\n', GammaII(i));
  end
  fprintf(fid, 'Eh = %.12f\n', Eh(i));
end
fprintf(fid, 'Boundary conditons\n');
for i=1:VesNum
  switch BCType(i,1)
    case 'u'
      fprintf(fid, 'u 0     # A lhs boundary Domain 1\n');
      fprintf(fid, '   a = %.12E\n', A(i));
      fprintf(fid, 'u 2     # U lhs boundary\n');
      fprintf(fid, '   u = %.12f\n', Vel(i));
      fprintf(fid, 'h 0     # Hd lhs boundary\n');
      fprintf(fid, '   h = 0.45\n');
    case 'p'
      fprintf(fid, 'p %d     # P lhs boundary Domain %d\n', Bifur(i,2), i);
      fprintf(fid, '   q = %.12f\n', 0);
      fprintf(fid, 'p 0     # U lhs boundary Domain\n');
      fprintf(fid, '   u = %.20E\n', 0);
      fprintf(fid, 'h 0     # Hd lhs boundary\n');
      fprintf(fid, '   h = %f\n', 0.45);
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
    case 'T'
      fprintf(fid, 'T  %.12f\n', 0.8);
      fprintf(fid, 'T  %.12f\n', 0.8);
      fprintf(fid, 'H  0\n');
    case 'u'
      fprintf(fid, 'u 0     # A lhs boundary Domain 84\n');
      fprintf(fid, '   a = %.12E\n', A(i));
      fprintf(fid, 'u 2     # U lhs boundary\n');
      fprintf(fid, '   u = %.12f\n', Vel(i));
      fprintf(fid, 'H 0     # Hd lhs boundary\n');
  end
end
fprintf(fid, 'Initial condition\n');
for i=1:VesNum
  fprintf(fid, 'a = Ao\n');
  fprintf(fid, 'u = %.12f\n', 0);
  fprintf(fid, 'h = %.12f\n', Hd(i));
end
fprintf(fid, 'History Pts\n');
fprintf(fid, '%d #Number of Domains with history points\n', VesNum);
for i=1:VesNum
  fprintf(fid, '3 %d  #Npts Domain id x[1], x[2], x[3]\n', i);
  fprintf(fid, '%.12f %.12f %.12f\n', Len(i)*0.0, Len(i)*0.5, Len(i)*1.0);
end
% for i=1:VesNum
%     fprintf(fid, '1 %d  #Npts Domain id x[1], x[2], x[3]\n', i);
%     fprintf(fid, '%.12f\n', Len(i)*0.5);
% end
fclose(fid);