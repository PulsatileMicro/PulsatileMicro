%% 生成输入文件
function fileName=GenInputFile_Systemic(VesType,VesParam,BCType,Bifur,BCVal,ModelParam)
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
TaperRate=VesParam(21,1);

dt=ModelParam(1);
Period=ModelParam(2);
NumCycles=ModelParam(3);
ODESolver=ModelParam(4);

VesNum=length(Len);
fileName=[VesType '.in'];
fid=fopen(fileName, 'w');

fprintf(fid, '15 parameter list\n');
fprintf(fid, '%E\tDT\n', dt);
fprintf(fid, '%d\tNSTEPS\n', uint32(Period/dt*NumCycles));
fprintf(fid, '%lu\tHISSTEP\n', uint32(1/dt/1000));
fprintf(fid, '%E\tRho (kg/m3)\n',1050);
fprintf(fid, '%f\tAlpha (velocity profile param in the momentum eq) (for Alpha=1 no wall viscous effects considered)\n', alpha(1));
fprintf(fid, '%f\tpinf(Pa)\n', 0);
fprintf(fid, '%d\tUpdVisc\n', 0);
fprintf(fid, '%d\tRiemann\n', 1);
fprintf(fid, '%d\tBinout\n', 0);
fprintf(fid, '%E\tRELTOL\n',1e-6);
fprintf(fid, '%E\tABSTOL\n',1e-6);
fprintf(fid, '%d\tSHOWCFL\n',0);
fprintf(fid, '%d\tSTEPLAPSE\n',1);
fprintf(fid, '%d\tODESolver\n',ODESolver);
fprintf(fid, '%d\tTAPER\n',0);

%%%% Mesh section %%%%
fprintf(fid, 'Mesh  -- expansion order -- quadrature order Ndomains = %d\n', VesNum);
for i=1:VesNum
  [A(i) Eh(i)]=Eval_Eh_A(Diam(i), E(i), WallTh(i));
%   if Gamma(1) ~= 0
%     fprintf(fid, '1\tnel Eh Area Visc Hd Gamma\n');
%   elseif GammaII(1) ~= 0
%     fprintf(fid, '1\tnel Eh Area Visc Hd GamII\n');
%   else
    fprintf(fid, '1\tnel Eh Area Visc Hd\n');
%   end
  fprintf(fid, '0.0\t %.12f \t%d\t%d\t# x_lower x_upper L q\n', Len(i), L, q);
  fprintf(fid, 'Ao = %.16E\n', A(i));
  fprintf(fid, 'Visc = %.12f\n', visc(i));
  fprintf(fid, 'Hd = %.12f\n', Hd(i));
  if Gamma(1) ~= 0
    fprintf(fid, 'Gamma = %f\n', Gamma(i));
  elseif GammaII(1) ~= 0
    fprintf(fid, 'GamII = %f\n', GammaII(i));
  end
  fprintf(fid, 'Eh = %.12f\n', Eh(i));
end
fprintf(fid, 'Boundary conditons\n');
for i=1:VesNum
  switch BCType(i,1)
    case 'q'
      fprintf(fid, 'q 0     # A lhs boundary Domain %d\n',i);
      fprintf(fid, '   a = %.20E\n', A(i));
      fprintf(fid, 'q %d     # Q lhs boundary\n', Bifur(i,2));
      fprintf(fid, '   q = %.20E\n', -1);  % Flow in mm^3/s
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
      fprintf(fid, 'R  %.16E\n', Bifur(i,3));
      fprintf(fid, 'R  %.16E\n', Bifur(i,3));
      fprintf(fid, 'H  0\n');
    case 'w'
      fprintf(fid, 'w %.16E\n', Bifur(i,3));
      fprintf(fid, 'w %.16E\n', Bifur(i,4));
      fprintf(fid, 'H  0\n');
    case 'W'
      fprintf(fid, 'W %.16E\n', Bifur(i,3));
      fprintf(fid, 'W %.16E\n', Bifur(i,4));
      fprintf(fid, 'H  0\n');
    case 'T'
      fprintf(fid, 'T  %f\n', Bifur(i,3));
      fprintf(fid, 'T  %f\n', Bifur(i,3));
      fprintf(fid, 'H  0\n');
    case 't'
      fprintf(fid, 't  %f\n', Bifur(i,3));
      fprintf(fid, 't  %f\n', Bifur(i,3));
      fprintf(fid, 'H  0\n');
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