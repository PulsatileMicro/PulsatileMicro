%% 生成模型流入边界所需的数据
% u - 流速，q - 流量，p - 压力
% BoundData格式: [边界值(u/q/p) 血管序号(1-546) Art/Ven标记(0-Art,1-Ven)]
function [i_artinput i_veninput]=...
  GenBoundInput(inFileName,outFileName,BoundData,ModelParam,VesParam,...
  nFreqTimes,PulsatLevel,dt,Period,PeriodNum,len_ratio,mass_ratio)
% 全局变量声明（定义在Macro.m中）
global NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE IN_BOUND_TYPE RANDOM_PHASE VEL_PROFILE USE_OPTBOUND BIFURLOSS
global HUMAN CATMEN HALFSIN IMPULSE INCRAMP DECRAMP SINE MULTI_FREQ_SIN SQUARE HUMAN_ROUTEA_546 HYPER VIT
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP
RandomPhaseFlag=ModelParam(RANDOM_PHASE);
Diam=VesParam(2,:);
scale_lamda=VesParam(17,1);
scale_u0=VesParam(18,1);
scale_r0=VesParam(19,1);
switch ModelParam(INPUT_WAVE)   % Input Waveform Type
  case HUMAN  % Sublingual
    % Column 1: Time (s)
    % Column 2: Velocity (um/s)
    % Column 3: Smoothed velocity (um/s)
    load 'HumanArt_Vel.mat';
    load 'HumanVen_Vel.mat';
    time=arteriole(:,1);
    artinput=arteriole(:,3)/1e3;  % (mm/s)
    veninput=venule(:,3)/1e3; % (mm/s)
    sel_ind=6:26;   % 选用的波形的下标，20 intervals, 20*0.04=0.8s/period
    % sel_ind=26:45;    % second pulse
  case CATMEN  % Cat mesenteric
    % delta_t: dt of the Gaehtgens data
    % mean_vel: mean velocity
    % art_vel: velocity series
    load 'gaehtgens_data.mat';
    time=0:delta_t:delta_t*26;time=time';
    artinput=art_vel/1e3;   % (mm/s)
    veninput=art_vel/1e3;
    sel_ind=1:27;
  case HALFSIN  % Half sine
    time=0:dt:Period;time=time';
    artinput=1/0.7*(1.575052015*sin(2*pi*time+0.1875*pi)-0.875052015)...
      .*((1.575052015*sin(2*pi*time+0.1875*pi)-0.875052015)>0);
    veninput=artinput;
    sel_ind=1:length(time);
  case IMPULSE  % Impulse
    % Refer to Alastruey's Doctoral Thesis
    time=0:dt:Period;time=time';
    t0=0.05;tau=0.002;
    artinput=1*exp(-((time-t0)/tau).^2);
    veninput=artinput;
    sel_ind=1:length(time);
  case INCRAMP  % Increasing Ramp
    time=0:dt:Period;time=time';
    artinput=zeros(length(time),1);
    ind=51:300;
    artinput(ind)=time(ind)-time(ind(1));
    veninput=artinput;
    sel_ind=1:length(time);
  case DECRAMP  % Decreasing Ramp
    time=0:dt:Period;time=time';
    artinput=zeros(length(time),1);
    ind=51:300;
    artinput(ind)=time(ind(end))-time(ind);
    veninput=artinput;
    sel_ind=1:length(time);
  case 6  % Sinusoidal
    time=0:dt:Period;time=time';
    artinput=sin(2*pi*time/Period*nFreqTimes)+2;
    veninput=artinput;
    sel_ind=1:length(time);
  case 7  % Multi-freq sine
    time=0:dt:Period;time=time';
    artinput=sin(2*pi*time/Period)+0.5*sin(2*pi*time/Period*2)+2;
    veninput=artinput;
    sel_ind=1:length(time);
  case 8  % Square
    time=0:dt:Period;time=time';
    artinput=square(2*pi*time/Period,30)+2;
    artinput(end)=artinput(end-1);
    veninput=artinput;
    sel_ind=1:length(time);
  case 9  % Vel on Route A in 546 network
    load 546_U_RouteA.mat;
    artinput=VelInput(PulsatLevel,:);
    artinput=[artinput artinput(1)]';
    veninput=VelInput(9,:);
    veninput=[veninput veninput(1)]';
    time=0:dt:Period*nFreqTimes;time=time';
    sel_ind=1:length(time);
  case 10  % Pressure from Ursino model
%     load NormPsp.mat;
    load PeriPsp.mat;
%     load AthePsp.mat;
%     NormPsp=AthePsp;
    NormPsp=PeriPsp;
    t1=0:1e-2:0.83;t2=0:1e-3:Period;
    artinput=interp1(t1,NormPsp,t2,'pchip')';
    veninput=artinput;
    time=0:dt:Period;time=time';
%     artinput=interp()
    sel_ind=1:length(time);
  case 11 % Egg818 Vel. from Exp.
    load Egg818_Vel.mat;
    artinput=arteriole;
    veninput=venule;
    time=time';
    sel_ind=1:length(time);
end

% one_xx: 一个周期的序列
% all_xx: N个周期的序列
one_time=time(sel_ind)-time(sel_ind(1));
one_time=one_time/nFreqTimes;
one_artinput=artinput(sel_ind);
one_veninput=veninput(sel_ind);
all_time=[];
all_artinput=[];
all_veninput=[];
    
for i=1:PeriodNum
  if i==1
    all_time=[all_time;one_time];
    all_artinput=[all_artinput;one_artinput];
    all_veninput=[all_veninput;one_veninput];
  else
    all_time=[all_time;one_time(2:end)+max(all_time)];
    all_artinput=[all_artinput;one_artinput(2:end)];
    all_veninput=[all_veninput;one_veninput(2:end)];
  end
end

[m,n]=size(BoundData);
for i=1:m
  VesNum=BoundData(i,2);
  % Scale the input waveform
  if BoundData(i,3)==0  % 动脉
    ratio=BoundData(i,1)/mean(all_artinput);
    all_artinput=all_artinput*ratio;
  else  % 静脉
    ratio=BoundData(i,1)/mean(all_veninput);
    all_veninput=all_veninput*ratio;
  end
  
  if ModelParam(INPUT_WAVE)==HUMAN || ModelParam(INPUT_WAVE)==CATMEN || ModelParam(INPUT_WAVE)==HUMAN_ROUTEA_546 ...
    || ModelParam(INPUT_WAVE)==HYPER || ModelParam(INPUT_WAVE)==VIT
    % 如果输入波形为实测波形，则进行插值，使采样频率变为1/dt Hz
    itime=0:dt:max(all_time);
    switch ModelParam(IN_BOUND_TYPE)
      case 0  % u
        i_artinput=interp1(all_time,all_artinput,itime,'pchip')/1e3/scale_u0;
        i_veninput=interp1(all_time,all_veninput,itime,'pchip')/1e3/scale_u0;
      case 1  % q
        if ModelParam(ODESOLVER)==RLC_EXP || ModelParam(ODESOLVER)==RLC_IMP ...
            || ModelParam(ODESOLVER)==RC_EXP || ModelParam(ODESOLVER)==RC_IMP
          % 如果是0D模型（RLC、RC），则需要进行量纲调整
          i_artinput=interp1(all_time,all_artinput,itime,'pchip')*pi*(Diam(VesNum)/2/len_ratio)^2*1e-3*1e-6*1e9*len_ratio^3;
          i_veninput=interp1(all_time,all_veninput,itime,'pchip')*pi*(Diam(VesNum)/2/len_ratio)^2*1e-3*1e-6*1e9*len_ratio^3;
        else
          i_artinput=interp1(all_time,all_artinput,itime,'pchip')/1e3*pi*(Diam(VesNum)/2)^2/scale_u0/(scale_r0^2);
          i_veninput=interp1(all_time,all_veninput,itime,'pchip')/1e3*pi*(Diam(VesNum)/2)^2/scale_u0/(scale_r0^2);
        end
      case 2  % p
        % TODO: 去量纲化
        i_artinput=interp1(all_time,all_artinput,itime,'pchip')*133;
        i_veninput=interp1(all_time,all_veninput,itime,'pchip')*133;
    end
  else
    % 如果输入波形为人工生成的波形，则不需要插值
    switch ModelParam(IN_BOUND_TYPE)
      case 0
        i_artinput=all_artinput/1e3/scale_u0;
        i_veninput=all_veninput/1e3/scale_u0;
      case 1
        i_artinput=all_artinput/1e3*pi*(Diam(VesNum)/2)^2/scale_u0/(scale_r0^2);
        i_veninput=all_veninput/1e3*pi*(Diam(VesNum)/2)^2/scale_u0/(scale_r0^2);
      case 2
        % TODO: 去量纲化
        i_artinput=all_artinput*133;
        i_veninput=all_veninput*133;
    end
  end
  
  % Add phase difference
  if i~=1 && RandomPhaseFlag==1
    i_artinput=circshift(i_artinput',round(100*rand))';
  end
  % Add frequent components
  % 用于测试粘弹性血管壁效应
%   f=32;
%   amp=mean(i_artinput);
%   y=amp/5*sin(2*pi*f*itime);
%   i_artinput=i_artinput+y;
  
  % 输出结果到.bcs文件
  if(BoundData(i,3)==0)  % inlet
    if length(Diam)==1
      fid=fopen([inFileName '.bcs'], 'w');
    else
      fid=fopen([inFileName '_' int2str(VesNum) '.bcs'], 'w');
    end
    fwrite(fid,i_artinput,'double',0,'l');
    fclose(fid);
    %         clear art_ivel;
  else    % outlet
    if length(Diam)==1
      fid=fopen([outFileName '.bcs'], 'w');
    else
      fid=fopen([outFileName '_' int2str(VesNum) '.bcs'], 'w');
    end
    fwrite(fid,i_veninput,'double',0,'l');
    fclose(fid);
  end
end
