%% 根据NetTypeID获取NetTypeName（字符串）
function NetTypeName=GetNetTypeName(NetTypeID)
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID SymNet_ID Junc_ID Single_ID Tree_ID
% NetTypeName用于给输出的文件命名
switch NetTypeID
  case Net_546_ID
    NetTypeName='Net_546';
  case Net_546_Meas_ID
    NetTypeName='Net_546_Meas';
  case Egg_818_ID
    NetTypeName='Egg_818';
  case Net_122_ID
    NetTypeName='Net_122';
  case Net_389_ID
    NetTypeName='Net_389';
  case Net_913_ID
    NetTypeName='Net_913';
  case Egg_CAM_ID
    NetTypeName='Egg_CAM';
  case Sub_CAM_ID
    NetTypeName='Sub_CAM';
  case Egg_636_ID
    NetTypeName='Egg_636';
  case SymNet_ID
    NetTypeName='SymNet';
  case Tree_ID
    NetTypeName='Tree';
  case Junc_ID
    NetTypeName='Junc';
  case Single_ID
    NetTypeName='Single';
end