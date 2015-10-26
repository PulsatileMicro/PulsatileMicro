function [NewAdapCoeff,ObjValue,OptDiam,OptWallTh,Sm,Sc] = STRUCT_ADAP_FUNC(NetTypeID,DataArray,Boundary,FuncPara,DatMatrix)
global WITH_WALL WITHOUT_WALL WITH_WALL_Cx WITHOUT_WALL_Cx WITH_WALL_Cx_PULSE
global PSO DOWNHILL GLOBAL_SEARCH YSPSO
global NOT_OPT_PARA OPT_PARA
global Net_546_ID Net_546_Meas_ID Egg_818_ID Net_122_ID Net_389_ID Net_913_ID Egg_CAM_ID Sub_CAM_ID Egg_636_ID

%% 1. ��������Ӧ�������Ż��������Ż�����
AdapType=WITH_WALL;
OptMethod=PSO;
OptType=OPT_PARA;

%% 2. ����Ԥ����
% 2.1 ���ݽ��
From=DatMatrix(:,2);
To=DatMatrix(:,3);
VesNum=length(From);

% 2.2 ��λ����
% NOTICE: �������ӳ���ͬ���˴���λ������mm�������������ĳ���
% Len=Len*1e-3;       % mm
% Diam=Diam*1e-3;     % mm
% Visc=Visc*1e-3;     % Pa.s
% WallTh=WallTh*1e-3; % mm

% 2.3 �߽緽������
% ���ݱ߽������ķ�������
% Ѫ��Ϊ����->��߽�->�߽�ڵ���From����
% Ѫ��Ϊ����->���߽�->�߽�ڵ���To����
[From,To]=OriNodeModify(Boundary,From,To);

%% 3. �������������ʼ����Ӧ
% 3.1 ����Ӧ��������
dt=0.1;       % ����Ӧ����
M0=1000;
Qref=0.001;   % ������������������ģ�飬�����Ż�
PO2ref=94.4;   % ����ѹ���ղ���������ģ��
J0=6618;       % �����źż������������ģ��
AdapPara=[dt,M0,Qref,PO2ref,J0];

switch AdapType
  case WITH_WALL    % 12������
    DatMatrix(:,6)=4e-3*ones(VesNum,1);
    load AdapCoeff_Wall.mat;
%     load AdapCoeff_Wall_Pries.mat;
    AdapCoeff=AdapCoeff_Wall;
  case WITHOUT_WALL
    % TODO
  case {WITH_WALL_Cx,WITH_WALL_Cx_PULSE}
    DatMatrix(:,6)=4e-3*ones(VesNum,1);
    load AdapCoeff_Wall.mat;
    load AdapCoeff_Wall_Cx.mat;
    AdapCoeff=AdapCoeff_Wall_Cx;
  case WITHOUT_WALL_Cx
    % TODO
end

if NetTypeID==Net_546_ID || NetTypeID==Net_122_ID
  % ��ȡFuncPara���������Ϊ�߽�
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,J0,FuncPara);
else
  % ��FuncPara
  [BHd,BSO2in,BJm,BJc]=BoundaryInput(Boundary,From,To,Qref,J0);
end
AdapBoundary=[BHd,BSO2in,BJm,BJc];

% ��ʼHd����˳��. Porder, ����. Norder, ����
Porder=1:VesNum;
Norder=VesNum:-1:1;
HdOrder=[Porder;Norder];

if OptType==OPT_PARA
  % ѡ���Ż��㷨
  switch OptMethod
    case PSO
      % TODO: �������������Ҫ����
      SwarmSize=10;
      InitSwarm=zeros(SwarmSize,length(AdapCoeff));
      
      %  'OutputFcns',@OptOutFunc,'PlotFcns',@OptPlotFunc
%         kc=AdapCoeff(1); 1.66��5
%         kmd=AdapCoeff(2); 0.955��5
%         kmg=AdapCoeff(3); -0.374��5
%         ksd=AdapCoeff(4); 3.077��5
%         ksg=AdapCoeff(5); 0.0177��5
%         kwt=AdapCoeff(6); 0.114��5
%         kwo=AdapCoeff(7); 0.609��5
%         tauref=AdapCoeff(8); 0.5598��80%
%         J0=AdapCoeff(9); 6618��20%
%         Lref=AdapCoeff(10); 14292��20%
%         Oref=AdapCoeff(11); 32050��80%
%         wref=AdapCoeff(12); 0.804��80%
        kcScale=linspace(1.66-5,1.66+5,SwarmSize);
        kmdScale=linspace(0.955-5,0.955+5,SwarmSize);
        kmgScale=linspace(-0.374-5,-0.374+5,SwarmSize);
        ksdScale=linspace(3.077-5,3.077+5,SwarmSize);
        ksgScale=linspace(0.0177-5,0.0177+5,SwarmSize);
        kwtScale=linspace(0.114-5,0.114+5,SwarmSize);
        kwoScale=linspace(0.609-5,0.609+5,SwarmSize);
        taurefScale=linspace(0.5598*0.2,0.5598*1.8,SwarmSize);
        J0Scale=linspace(6618*0.8,6618*1.2,SwarmSize);
        LrefScale=linspace(14292*0.8,14292*1.2,SwarmSize);
        OrefScale=linspace(32050*0.2,32050*1.8,SwarmSize);
        wrefScale=linspace(0.804*0.2,0.804*1.8,SwarmSize);
      for i=1:SwarmSize
 %       InitSwarm(i,:)=AdapCoeff.*(rand(1,length(AdapCoeff))+0.5);
        InitSwarm(i,:)=[kcScale(1,i),kmdScale(1,i),kmgScale(1,i),ksdScale(1,i),ksgScale(1,i),kwtScale(1,i),kwoScale(1,i),taurefScale(1,i),J0Scale(1,i),LrefScale(1,i),OrefScale(1,i),wrefScale(1,i)];
      end
      lb=[1.66-5,0.955-5,-0.374-5,3.077-5,0.0177-5,0.114-5,0.609-5,0.5598*0.2,6618*0.8,14292*0.8,32050*0.2,0.804*0.2];
      ub=[1.66+5,0.955+5,-0.374+5,3.077+5,0.0177+5,0.114+5,0.609+5,0.5598*1.8,6618*1.2,14292*1.2,32050*1.8,0.804*1.8];
      options = optimoptions('particleswarm','SwarmSize',SwarmSize,'InitialSwarm',InitSwarm,'MaxIter',1000,'UseParallel',false,'OutputFcns',@pswIteraRes);
      [X,FVAL,EXITFLAG,OUTPUT]= particleswarm(@(x) AdapObjFunc(x,NetTypeID,AdapType,HdOrder,...
        AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray),length(AdapCoeff),lb,ub,options);
      NewAdapCoeff=X;
%   EXITFLAG 
%   1   Relative change in the objective value over the last options.StallIterLimit iterations is less than options.TolFun.
%   0   Number of iterations exceeded options.MaxIter.
%   -1  Iterations stopped by output function or plot function.
%   -2  Bounds are inconsistent: for some i, lb(i) > ub(i).
%   -3  Best objective function value is at or below options.ObjectiveLimit.
%   -4  Best objective function value did not change within options.StallTimeLimit seconds.
%   -5  Run time exceeded options.MaxTime seconds.
      formatstring = 'particleswarm reached the value %f using %d function evaluations, after %d iterations with ExitFlag %d and algorithm stop reason is %s.\n';
      fprintf(formatstring, FVAL, OUTPUT.funccount, OUTPUT.iterations, EXITFLAG, OUTPUT.message);
    case YSPSO
        SwarmSize = 1;
        InitSwarm=zeros(SwarmSize,length(AdapCoeff));
        for i=1:SwarmSize
          InitSwarm(i,:)=AdapCoeff.*(rand(1,length(AdapCoeff))+0.5);
        end
        c1 = 2.8;
        c2 = 1.3;
        AdaptTimes = 1;
        [X,FVAL] = YasuoPSO(@(x) AdapObjFunc(x,NetTypeID,AdapType,HdOrder,...
          AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray),SwarmSize,InitSwarm,c1,c2,AdaptTimes,length(AdapCoeff));
      NewAdapCoeff=X;
    case DOWNHILL
      % Simplex Downhill����
      options=optimset('TolFun',1e-2,'TolX',1e-2,'MaxIter',100,'Display','iter');
      [X,FVAL,EXITFLAG,OUTPUT]= fminsearch(@(x) AdapObjFunc(x,NetTypeID,AdapType,HdOrder,...
        AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray),AdapCoeff,options);
      NewAdapCoeff=X;
%   EXITFLAG 
%   1   Relative change in the objective value over the last options.StallIterLimit iterations is less than options.TolFun.
%   0   Number of iterations exceeded options.MaxIter.
%   -1  Iterations stopped by output function or plot function.
%   -2  Bounds are inconsistent: for some i, lb(i) > ub(i).
%   -3  Best objective function value is at or below options.ObjectiveLimit.
%   -4  Best objective function value did not change within options.StallTimeLimit seconds.
%   -5  Run time exceeded options.MaxTime seconds.
      formatstring = 'Nelder-Mead reached the value %f using %d function evaluations, after %d iterations with ExitFlag %d and algorithm stop reason is %s.\n';
      fprintf(formatstring, FVAL, OUTPUT.funcCount, OUTPUT.iterations, EXITFLAG, OUTPUT.message);
    case GLOBAL_SEARCH
      % TODO
  end
  % ʹ�ø��º��ϵ����һ��
  [ObjValue,AdapCoeff,OptDiam,OptWallTh,Sm,Sc]=AdapObjFunc(NewAdapCoeff,NetTypeID,AdapType,HdOrder,...
    AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray);
  NewAdapCoeff=AdapCoeff;
else
  [ObjValue,AdapCoeff,OptDiam,OptWallTh,Sm,Sc]=AdapObjFunc(AdapCoeff,NetTypeID,AdapType,HdOrder,...
    AdapBoundary,AdapPara,Boundary,DatMatrix,DataArray);
  NewAdapCoeff=AdapCoeff;
end