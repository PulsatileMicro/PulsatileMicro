% 剪切力信号、压力信号 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Stau,Sp,Tau]=SCounter(DeltaP,Diam,Len,MeanP,Tauref,Eju)
% shear stress:t=D*P/(4L);
if Eju==0  %非正常情况
    % Stau
    Tau=133.*abs(DeltaP).*1e1.*Diam./(4.*Len);
    Stau=log10(abs(Tau)+Tauref);
    % Sp
%     InP=MeanP./1.332;   %mmHg 不知为何，需要除以1.33才能得到与数据一致的结果
    InP=MeanP; 
    InPIndex= InP<10;   %异常处理：小于10mmHg，Sp计算会出错（复数）
    InP(InPIndex)=10;   %Pressure<10mmHg,则Tau为14dyn/cm2
    TauP=100-86.*exp(-5000.*((log10(log10(InP))).^5.4));
    Sp=-log10(TauP);
else
    VelNum=length(Diam);
    Tau=zeros(VelNum,1);
    Stau=zeros(VelNum,1);
    Sp=zeros(VelNum,1); 
end