%% ���ݷ������ȷ�������ź�����
function Period=GetPeriod(NetTypeID,ModelParam)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID SymNet_ID Junc_ID Single_ID Tree_ID
global NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE IN_BOUND_TYPE RANDOM_PHASE VEL_PROFILE USE_OPTBOUND BIFURLOSS
global HUMAN CATMEN HALFSIN IMPULSE INCRAMP DECRAMP SINE MULTI_FREQ_SIN SQUARE HUMAN_ROUTEA_546 HYPER VIT

% ��ͬ���ε��Ķ�����ʱ�䣬��λ����
CAT_MEN_PERIOD=0.325;               % Gaehtgens��õ�è��ϵĤ���ٲ��Σ�����Ϊ0.325��
VIT_PERIOD=0.28;                    % �������벨������
HUMAN_PERIOD=0.8;

% è��ϵĤ��������(ͨ�����ڴ���ϵĤ)
if ModelParam(INPUT_WAVE)==CATMEN
  Period=CAT_MEN_PERIOD;
% ���߲�������
elseif NetTypeID==Egg_818_ID || NetTypeID==Egg_CAM_ID ||...
    NetTypeID==Sub_CAM_ID || NetTypeID==Egg_636_ID
  Period=VIT_PERIOD;
  ModelParam(INPUT_WAVE)=VIT;        
elseif ModelParam(INPUT_WAVE)==SINE   % ����������ź�����
  Period=1;                           % �����ź�����ʱ�����ڣ���������1,2,4,8,16,32Hz������
else
  Period=HUMAN_PERIOD;                % �������
end
