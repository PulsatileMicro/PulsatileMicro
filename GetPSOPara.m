function [psoparams,PSOseedValue] = GetPSOPara(AdapCoeff)
shw     = 0;
ps      = 2;
ac      = [2.1,2.1];% acceleration constants, only used for modl=0
Iwt     = [0.9,0.6];  % intertia weights, only used for modl=0
epoch   = 400; % max iterations
wt_end  = 100; % iterations it takes to go from Iwt(1) to Iwt(2), only for modl=0
errgrad = 1e-99;   % lowest error gradient tolerance
errgraditer=100; % max # of epochs without error change >= errgrad
errgoal = NaN;
%                 0 = Common PSO w/intertia (default)
%                 1,2 = Trelea types 1,2
%                 3   = Clerc's Constricted PSO, Type 1"
modl    = 0;    
PSOseed = 1;    % if=1 then can input particle starting positions, if= 0 then all random

PSOseedValue=zeros(ps,length(AdapCoeff));
% 初始化PSO种子
for i=1:ps
%   PSOseedValue(i,:)=AdapCoeff*5.*rand(1,length(AdapCoeff));
  PSOseedValue(i,:)=AdapCoeff;
end

psoparams=...
  [shw epoch ps ac(1) ac(2) Iwt(1) Iwt(2) ...
  wt_end errgrad errgraditer errgoal modl PSOseed];
end
