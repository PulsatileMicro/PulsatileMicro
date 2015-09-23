% Read in the velocity data and generate the input file for 1D model
function []=GaussVelProc(dt,inFileName,outFileName,Vel,periodNum,scale_u0)

% Column 1: Time
% Column 2: Velocity
% Column 3: Smoothed velocity

% dt=dt*1000;
t=0:dt:0.8*(periodNum+1);
t0=0.05;
tau=0.01;

[m,n]=size(Vel);
for i=1:m
  VesNum=Vel(i,2);
  % Normalize the velocity to the value from the data file
%   y=t.*(t<0.3);
  if Vel(m,3)==0
    art_ivel=Vel(i,1)*exp(-((t-t0)/tau).^2)/1e6;
%     art_ivel=Vel(i,1).*y/1e6;
  else
    
  end
  
  if Vel(i,3)==0  % inlet
    if VesNum==1
      fid=fopen(inFileName, 'w');
    else
      fid=fopen([inFileName int2str(VesNum) '.bcs'], 'w');
    end
    fwrite(fid,art_ivel,'double',0,'l');
    fclose(fid);
  else
  end
end