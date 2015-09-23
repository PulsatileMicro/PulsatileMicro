%% Adjust Diameter
% ���ڹܾ������ڲ��Թܾ���С��Damping��Ӱ��
% �Ŵ�ܾ�ʱ������ëϸѪ�ܳߴ粻��
% ��С�ܾ�ʱ���������΢�����ߴ粻��
function [Diam DiamRatio]=AdjustDiam(Diam,VesType,AdjMode)
% AdjMode �ܾ�����ģʽ
% 0: ���ֲ���
% 1�����ԷŴ�2��ָ���Ŵ�3�������Ŵ�
% 4��������С��5��ָ����С��6��������С
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
