% Read in the velocity data and generate the input file for 1D model
function MeanBottomRatioA=VelProc(dt,inFileName,outFileName,artVel,venVel,periodNum,scale_u0)
% Column 1: Time
% Column 2: Velocity
% Column 3: Smoothed velocity
load 'arteriole.mat';
load 'venule.mat';
time=arteriole(:,1);
art_vel=arteriole(:,3);
ven_vel=venule(:,3);

t_one_pulse=6:26;   % 20 intervals, 20*0.04=0.8s/period first pulse
% t_one_pulse=26:45;    % second pulse
new_time=time(t_one_pulse)-time(t_one_pulse(1));
new_art_vel=art_vel(t_one_pulse);
new_ven_vel=ven_vel(t_one_pulse);
new_time_all=[];
new_artvel_all=[];
new_venvel_all=[];

for i=1:periodNum+1 % Don't know why!!
    if i==1
        new_time_all=[new_time_all;new_time];
        new_artvel_all=[new_artvel_all;new_art_vel];
        new_venvel_all=[new_venvel_all;new_ven_vel];
    else
        new_time_all=[new_time_all;new_time(2:end)+max(new_time_all)];
        new_artvel_all=[new_artvel_all;new_art_vel(2:end)];
        new_venvel_all=[new_venvel_all;new_ven_vel(2:end)];
    end
    
end

% Normalize the velocity to the value from the data file
ratio=artVel/mean(new_artvel_all)*133;
new_artvel_all=new_artvel_all*ratio;
ratio=venVel/mean(new_venvel_all)*133;
new_venvel_all=new_venvel_all*ratio;
% The ratio between the mean and bottom value of the waveform
% Used for changing the mesentery diam (mean) to bottom value as the init value
MeanBottomRatioA=mean(new_artvel_all)/min(new_artvel_all);

% interpolation
% dt=dt*1000;
itime=0:dt:max(new_time_all);
art_ivel=interp1(new_time_all, new_artvel_all, itime, 'cubic')/scale_u0;
ven_ivel=interp1(new_time_all, new_venvel_all, itime, 'cubic')/scale_u0;

% Add frequent components
% f=32;
% amp=mean(art_ivel);
% y=amp/5*sin(2*pi*f*itime);
% art_ivel=art_ivel+y;
% clear itime;

% write solution
fid=fopen(inFileName, 'w');
fwrite(fid,art_ivel,'double',0,'l');
fclose(fid);
clear art_ivel;
fid=fopen(outFileName, 'w');
fwrite(fid,ven_ivel,'double',0,'l');
fclose(fid);
