%% 根据仿真参数确定输入信号周期
function Period=GetPeriod(NetTypeID,ModelParam)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID SymNet_ID Junc_ID Single_ID Tree_ID
global NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE IN_BOUND_TYPE RANDOM_PHASE VEL_PROFILE USE_OPTBOUND BIFURLOSS
global HUMAN CATMEN HALFSIN IMPULSE INCRAMP DECRAMP SINE MULTI_FREQ_SIN SQUARE HUMAN_ROUTEA_546 HYPER VIT

% 不同波形的心动周期时间，单位：秒
CAT_MEN_PERIOD=0.325;               % Gaehtgens测得的猫肠系膜流速波形，周期为0.325秒
VIT_PERIOD=0.28;                    % 鸡胚输入波形周期
HUMAN_PERIOD=0.8;

% 猫肠系膜波形输入(通常用于大鼠肠系膜)
if ModelParam(INPUT_WAVE)==CATMEN
  Period=CAT_MEN_PERIOD;
% 鸡胚波形输入
elseif NetTypeID==Egg_818_ID || NetTypeID==Egg_CAM_ID ||...
    NetTypeID==Sub_CAM_ID || NetTypeID==Egg_636_ID
  Period=VIT_PERIOD;
  ModelParam(INPUT_WAVE)=VIT;        
elseif ModelParam(INPUT_WAVE)==SINE   % 如果是正弦信号输入
  Period=1;                           % 正弦信号输入时的周期，便于设置1,2,4,8,16,32Hz的输入
else
  Period=HUMAN_PERIOD;                % 其它情况
end
