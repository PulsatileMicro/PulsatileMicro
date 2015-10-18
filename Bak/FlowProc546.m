% Read in the velocity data and generate the input file for 1D model
function MeanBottomRatioA=FlowProc546(dt,inFileName,outFileName,Vel,periodNum,Diam,scale_u0,scale_r0)

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
% % scale the time
% new_time=new_time/4;
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

[m,n]=size(Vel);
for i=1:m
    VesNum=Vel(i,2);
    % Normalize the velocity to the value from the data file
    if Vel(m,3)==0
        ratio=Vel(i,1)/mean(new_artvel_all)*1000;
        new_artvel_all=new_artvel_all*ratio;
    else
        ratio=Vel(i,1)/mean(new_venvel_all)*1000;
        new_venvel_all=new_venvel_all*ratio;
    end
    
    % interpolation
%     dt=dt*1000;
    itime=0:dt:max(new_time_all);
    art_ivel=interp1(new_time_all, new_artvel_all, itime, 'cubic')/1e6*pi*(Diam(VesNum)*1e-6).^2/4/scale_u0/scale_r0/scale_r0;
    ven_ivel=interp1(new_time_all, new_venvel_all, itime, 'cubic')/1e6*pi*(Diam(VesNum)*1e-6).^2/4/scale_u0/scale_r0/scale_r0;
    
    % Add frequent components
    % f=32;
    % amp=mean(art_ivel);
    % y=amp/5*sin(2*pi*f*itime);
    % art_ivel=art_ivel+y;
    clear itime;
    
    % write solution
    
    if(Vel(i,3)==0)  % inlet
        fid=fopen([inFileName int2str(VesNum) '.bcs'], 'w');
        fwrite(fid,art_ivel,'double',0,'l');
        fclose(fid);
%         clear art_ivel;
    else    % outlet
        fid=fopen([outFileName int2str(VesNum) '.bcs'], 'w');
%         fwrite(fid,ven_ivel,'double',0,'l');
        fwrite(fid,art_ivel,'double',0,'l');
        fclose(fid);
    end
end

% The ratio between the mean and bottom value of the waveform
% Used for changing the mesentery diam (mean) to bottom value as the init value
MeanBottomRatioA=mean(new_artvel_all)/min(new_artvel_all);