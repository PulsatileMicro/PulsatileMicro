% Read in the velocity data and generate the input file for 1D model
function []=VelProc(dt,inFileName,outFileName,Vel,periodNum,scale_u0)

% Column 1: Time
% Column 2: Velocity
% Column 3: Smoothed velocity

% dt=dt*1000;
t=0:dt:0.8*(periodNum+1);

[m,n]=size(Vel);
for i=1:m
  VesNum=Vel(i,2);
  % Normalize the velocity to the value from the data file
  if Vel(m,3)==0
    art_ivel = Vel(i,1)/1e6/0.7*(1.575052015*sin(2*pi*t/0.8+0.1875*pi)-0.875052015)...
      .*((1.575052015*sin(2*pi*t/0.8+0.1875*pi)-0.875052015)>0)/scale_u0;% write solution
    art_ivel=art_ivel+0.07*mean(art_ivel);
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