function FQe=PS_effect(Df,Da,Db,Hd,FQb)
%Df����Ѫ�ܹܾ���um
%Da����ǰ��Ѫ�ܹܾ���um
%Db����һ��Ѫ�ܹܾ���um
%Hd����Ѫ�ܺ�ϸ������
%FQb����ǰ��Ѫ�ܵ�Ѫ����������ǰ��Ѫ��Ѫ��/��Ѫ��Ѫ����
A=-13.29*((Da*Da-Db*Db)/(Da*Da+Db*Db))*(1-Hd)/Df;
B=1+6.98*(1-Hd)/Df;
X0=0.964*(1-Hd)/Df;
% if FQb<X0
%   FQe=0;
% elseif FQb>1-X0
%   FQe=1;
% else
%   FQe=1/(1+exp(-(A+B*logit((FQb-X0)/(1-2*X0)))));
% end
FQe=1/(1+exp(-(A+B*logit((FQb-X0)/(1-2*X0)))));

function y=logit(x)
y=log(x/(1-x));
