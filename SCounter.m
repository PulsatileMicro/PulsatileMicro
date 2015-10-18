% �������źš�ѹ���ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Stau,Sp,Tau]=SCounter(DeltaP,Diam,Len,MeanP,Tauref,Eju)
% shear stress:t=D*P/(4L);
if Eju==0  %���������
    % Stau
    Tau=133.*abs(DeltaP).*1e1.*Diam./(4.*Len);
    Stau=log10(abs(Tau)+Tauref);
    % Sp
%     InP=MeanP./1.332;   %mmHg ��֪Ϊ�Σ���Ҫ����1.33���ܵõ�������һ�µĽ��
    InP=MeanP; 
    InPIndex= InP<10;   %�쳣����С��10mmHg��Sp��������������
    InP(InPIndex)=10;   %Pressure<10mmHg,��TauΪ14dyn/cm2
    TauP=100-86.*exp(-5000.*((log10(log10(InP))).^5.4));
    Sp=-log10(TauP);
else
    VelNum=length(Diam);
    Tau=zeros(VelNum,1);
    Stau=zeros(VelNum,1);
    Sp=zeros(VelNum,1); 
end