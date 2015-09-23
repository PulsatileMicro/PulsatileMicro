function [t P U Q A Visc VesID]=readHisFile(fileName,NumHisPt,ODESolver)
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP
fid=fopen(fileName, 'r');
if ODESolver==ONED_EXP || ODESolver==ONED_IMP
  C=textscan(fid, '%f %f %f %f %f %f %d', 'HeaderLines',3+NumHisPt);
else
  C=textscan(fid, '%f %f %f %f %f %f %d', 'HeaderLines',2);
end
% C=textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %d', 'HeaderLines',3+NumHisPt);
fclose(fid);
t=C{1};
P=C{2};
U=C{3};
Q=C{4};
A=C{5};
Visc=C{6};
VesID=C{7};
end