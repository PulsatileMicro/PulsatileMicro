%%%% 仿真参数设置 %%%%
global NONDIM VISC_UPD RIEM BINOUT VISCOELAS CFL STEP_LAPSE ODESOLVER INPUT_WAVE IN_BOUND_TYPE RANDOM_PHASE LIN_PA VEL_PROFILE USE_OPTBOUND BIFURLOSS
NONDIM=2;           % Non-dimensionalization flag 1: Nondimensionalization
VISC_UPD=3;         % Viscosity update strategy 0: No viscosity update 1: vitro 2: ESL vivo 3: No ESL vivo
RIEM=4;             % Riemann solver flag 0: Not use 1: Use riemann solver for the BioFlux
BINOUT=5;           % Binary output flag 0: Output text format 1: Output binary format
VISCOELAS=6;        % Viscoelastic wall 0: No viscous wall 1: New viscoelastic model 2: Old viscoelastic model
CFL=7;              % Show CFL flag 0: Hide CFL number 1: Show CFL number
STEP_LAPSE=8;       % Show Step Lapse flag 0: Hide Step Lapse 1: Show Step Lapse
ODESOLVER=9;        % ODE solver flag 0: 1D_Exp 1: 1D_Imp 2: RLC_Exp 3: RLC_Imp 4: Steady State
                    % 5,6: Womersley (Not Ready) 7: RC_Exp 8: RC_Imp 
                    % 9: Sparse Steady State Solver
INPUT_WAVE=10;      % Input Waveform 0: Sublingual 1: Cat mesentery 2: Half sine 3:Impulse 
                    % 4: Increasing Ramp 5: Decreasing Ramp 6: Sinusoidal 7: Multi-freq sine 
                    % 8: Square 9: Sublingual Vel on Route A 10: Hypertension 0D-1D Coupling Input 
                    % 11: Egg818 measured velocity
IN_BOUND_TYPE=16;   % BoundType 0:u(mm/s) 1:q(ml/s) 2:p(mmHg)
RANDOM_PHASE=17;    % Introduce random phase difference to secondary boundaries (in 546 network): 0: No 1: Yes
LIN_PA=18;          % Linearity of Pressure-Area relationship: 0: Linear 1: Non-Linear
VEL_PROFILE=19;     % Velocity profile: 0: Parabolic 1: Measured
USE_OPTBOUND=22;    % Use optimized boundary 0:not use 1:use
BIFURLOSS=23;       % Energy loss at the bifurcations 0: no loss 1:loss

%%%% ODE Solver
% 0: 1D_Exp 1: 1D_Imp 2: RLC_Exp 3: RLC_Imp 4: Steady State
% 5,6: Womersley (Not Ready) 7: RC_Exp 8: RC_Imp
global ONED_EXP ONED_IMP RLC_EXP RLC_IMP SS WOM_1 WOM_2 RC_EXP RC_IMP Sparse_SS
ONED_EXP=0;
ONED_IMP=1;
RLC_EXP=2;
RLC_IMP=3;
SS=4;
WOM_1=5;
WOM_2=6;
RC_EXP=7;
RC_IMP=8;
Sparse_SS=9;

%%%% 输入边界类型 %%%%
% Input Waveform 0: Sublingual 1: Cat mesentery 2: Half sine 3:Impulse 
% 4: Increasing Ramp 5: Decreasing Ramp 6: Sinusoidal 7: Multi-freq sine 
% 8: Square 9: Sublingual Vel on Route A 10: Hypertension 0D-1D Coupling Input 
% 11: Egg818 measured velocity
global HUMAN CATMEN HALFSIN IMPULSE INCRAMP DECRAMP SINE MULTI_FREQ_SIN SQUARE HUMAN_ROUTEA_546 HYPER VIT
HUMAN=0;
CATMEN=1;
HALFSIN=2;
IMPULSE=3;
INCRAMP=4;
DECRAMP=5;
SINE=6;
MULTI_FREQ_SIN=7;
SQUARE=8;
HUMAN_ROUTEA_546=9;   % 用HUMAN波形作为输入进行仿真，取出通路A上典型血管段的波形作为输入。需load 546_U_RouteA.mat
HYPER=10;
VIT=11;

%%%% 血管网络类型 %%%%
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID SymNet_ID Junc_ID Single_ID Tree_ID Heart
Net_546_ID=1;       % 1 - 546肠系膜血管网络（自适应后）: Net_546
Net_546_Meas_ID=2;  % 2 - 546肠系膜血管网络（自适应前）：Net_546_Meas
Egg_818_ID=3;       % 3 - 818鸡胚血管网络: Egg_818
Net_122_ID=4;       % 4 - 122人工血管网络: Net_122
Net_389_ID=5;       % 5 - 389肠系膜网络：Net_389
Net_913_ID=6;       % 6 - 913肠系膜网络：Net_913
Egg_CAM_ID=7;       % 7 - CAM鸡胚血管网络：Egg_CAM
Sub_CAM_ID=8;       % 8 - CAM鸡胚血管网络子网络：Sub_CAM
Egg_636_ID=9;       % 9 - 636鸡胚血管网络：Egg_636
Heart=10;           % 10 - Coronary tree
SymNet_ID=100;      % 100 - 对称血管网络: SymNet
Tree_ID=101;        % 101 - 对称血管树: Tree
Junc_ID=102;        % 102 - 接合血管: Junc
Single_ID=103;      % 103 - 单段血管: Single

%%%% Damping因素排除模式: DampFactor
global DampInit E10 Visc01 E01 Visc10 EA10 EC10 EV10 EA01 EC01 EV01 ViscA10 ViscC10 ViscV10 ViscA01 ViscC01 ViscV01 ViscBD01 ViscBD10 DiamN DiamAN DiamCN DiamVN E5 E02 Visc5 Visc02 DampE DampVisc
DampInit=1;     % 1 - 初始模式
E10=2;          % 2 - 杨氏模量增大10倍
Visc01=3;       % 3 - 粘滞度减小10倍
E01=4;          % 4 - 杨氏模量减小10倍
Visc10=5;       % 5 - 粘滞度增大10倍
EA10=6;         % 6 - 微动脉杨氏模量增大10倍
EC10=7;         % 7 - 毛细血管杨氏模量增大10倍
EV10=8;         % 8 - 微静脉杨氏模量增大10倍
EA01=9;         % 9 - 微动脉杨氏模量减小10倍
EC01=10;        % 10 - 毛细血管杨氏模量减小10倍
EV01=11;        % 11 - 微静脉杨氏模量减小10倍
ViscA10=12;     % 12 - 微动脉粘滞度增大10倍
ViscC10=13;     % 13 - 毛细血管粘滞度增大10倍
ViscV10=14;     % 14 - 微静脉粘滞度增大10倍
ViscA01=15;     % 15 - 微动脉粘滞度减小10倍
ViscC01=16;     % 16 - 毛细血管粘滞度减小10倍
ViscV01=17;     % 17 - 微静脉粘滞度减小10倍
ViscBD01=18;    % 18 - 边界阻力减小10倍
ViscBD10=19;    % 19 - 边界阻力增大10倍
DiamN=20;       % 20 - 管径增大N倍
DiamAN=21;      % 21 - 微动脉管径增大1倍
DiamCN=22;      % 22 - 毛细血管管径增大1倍
DiamVN=23;      % 23 - 微静脉管径增大1倍
E5=24;          % 24 - 杨氏模量增大5倍
E02=25;         % 25 - 杨氏模量减小5倍
Visc5=26;       % 26 - 粘滞度增大5倍
Visc02=27;      % 27 - 粘滞度减小5倍
DampE=28;
DampVisc=29;

%%%% main程序运行模式: RunType
global SIM_PREP FREQ_ANAL BOUND_OPT STRUCT_ADAP
SIM_PREP=1;
FREQ_ANAL=2;
BOUND_OPT=3;
STRUCT_ADAP=4;

%%%% 结构自适应仿真模式：AdapType
global WITH_WALL WITHOUT_WALL WITH_WALL_Cx WITHOUT_WALL_Cx WITH_WALL_Cx_PULSE
WITH_WALL=1;          % 包含血管壁优化
WITHOUT_WALL=2;       % 不包含血管壁优化
WITH_WALL_Cx=3;       % 包含血管壁优化，优化Cx参数
WITHOUT_WALL_Cx=4;    % 不包含血管壁优化，优化Cx参数
WITH_WALL_Cx_PULSE=5; % 包含血管壁优化，优化Cx参数，脉动情况
%%%% 结构自适应优化分类: OptCategory
global PSO SIMPLEX GA
PSO=1;
SIMPLEX=2;
GA=3;
%%%% 结构自适应优化方法: OptMethod
global PSORESV DOWNHILL GLOBAL_SEARCH YSPSO SELPSO
PSORESV=1;              % Matlab自带粒子群优化方法
DOWNHILL=2;         % Simplex Downhill方法(matlab fminsearch函数)
GLOBAL_SEARCH=3;
YSPSO=4;
SELPSO=5;
%%%% 优化类型：是否估算系数？:OptType
global NOT_OPT_PARA OPT_PARA OPT_BOUND_FLOW OPT_BOUND_META
NOT_OPT_PARA=1;
OPT_PARA=2;
OPT_BOUND_FLOW=3;
OPT_BOUND_META=4;