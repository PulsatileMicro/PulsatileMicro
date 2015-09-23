%% Adjust Diameter
% 调节管径，用于测试管径大小对Damping的影响
% 放大管径时，保持毛细血管尺寸不变
% 缩小管径时，保持入口微动脉尺寸不变
function [Diam DiamRatio]=AdjustDiam(Diam,VesType,AdjMode)
% AdjMode 管径调节模式
% 0: 保持不变
% 1：线性放大，2：指数放大，3：对数放大
% 4：线性缩小，5：指数缩小，6：对数缩小
switch AdjMode
  case 0
    DiamRatio=ones(length(Diam),1);
  case 1
    DiamRatio=Diam./min(Diam);
  case 2
    DiamRatio=Diam./min(Diam);
    DiamRatio=1.4.^(DiamRatio-1);
  case 3
    DiamRatio=Diam./min(Diam);
    DiamRatio=logn(DiamRatio+0.15,1.15);
  case 4
    DiamRatio=Diam./max(Diam(VesType==1));
  case 5
    DiamRatio=Diam./max(Diam(VesType==1));
    DiamRatio=1.40.^(DiamRatio-1);
  case 6
    DiamRatio=Diam./max(Diam(VesType==1));
    DiamRatio=logn(DiamRatio+0.15,1.15);
end
Diam=Diam.*DiamRatio;
