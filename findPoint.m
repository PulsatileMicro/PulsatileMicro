function [FE,FP,PPG_C2p,PPG_C3p,PPG_Ap,PPG_A2p,PPG_A3p,PPG_A4p,PPG_Bp,Tbp,Tep,Tpeakp,Qpeakp,Rp,HR,PointDelay]=findPoint(t,E,P,nFreq,SP,b_for)
%first find the point of ECG
fcutoff1=12/nFreq/2;
fcutoff2=0.5/nFreq/2;
order_lp=2;
order_hp=2;
% [b,a]=butter(1,fcutoff1);  %12Hz��ͨ
% [d,c]=butter(1,fcutoff2,'high'); %��ͨ��1���ȣ�ȥ����Ư��
b_lp=fir1(order_lp,fcutoff1);
b_hp=fir1(order_hp,fcutoff2,'high');
FE=filter(b_lp,1,E);
FP=filter(b_lp,1,P);
FE=filter(b_hp,1,FE);
FP=filter(b_hp,1,FP);
PointDelay=fix((order_lp+order_hp)/2);

FE(1:2*PointDelay)=[];
FP(1:2*PointDelay)=[];
diffE=[0,diff(FE)];
diffE=filter(b_lp,1,diffE);
diffE=filter(b_hp,1,diffE);
diffE(1:PointDelay)=[];

diffFP=[0,diff(FP)]; % one order difference
diffFP=filter(b_lp,1,diffFP);
diffFP=filter(b_hp,1,diffFP);
diffFP(1:PointDelay)=[];
 
L= length(diffE);
FE=FE(1:L);
FP=FP(1:L);
t=t*nFreq;
t=t-t(1);
newT=t(1:L);
% Em=(FE-mean(FE))./std(FE); %��һ������
% [b,a]=butter(1,fcutoff1);  %12Hz��ͨ
% b=fir1(4,fcutoff1);
figure;
subplot(2,1,1)
plot(newT,FE)
grid on
hold on
count=0;
flag=0;
j=1;
[Rp,HR,Rdirection]=GetQRS(FE,nFreq);   %�ҵ����е�R�����λ��
Swin=floor(60*nFreq/HR);
Swin2=floor(Swin/2);
denTHQb=1.8;    %������ֵ��ĸ���趨
 while j<=length(Rp) & flag==0 
    %RxxΪ��ʱ�����ϵ�ʵ��λ��
    %Ѱ��dermax(����źŵ����ֵ��)
    %dermaxx Ϊdemax������λ��
     if Rp(j)>Swin2
         i=Rp(j)-Swin2:Rp(j)+Swin2;
     else
         i=1:Rp(j)+Swin2;
     end
     if i(length(i))>=L
         flag=1;
         break;
     else
         [dermay(j),dermaxx(j)]=max((diffE(i(1):i(length(i)))));
         dermaxx(j)=dermaxx(j)+i(1);
         dermay(j)=diffE(dermaxx(j));
         [derminy(j),derminx(j)]=min((diffE(i(1):i(length(i)))));
         derminx(j)=derminx(j)+i(1);
         derminy(j)=diffE(derminx(j));
     end
     count=count+1;
    %%find the S wave position
     pk_pre=fix(nFreq*20/1000);
     pk_back=fix(nFreq*120/1000);
     i=dermaxx(j)+pk_pre:dermaxx(j)+pk_back; %��λpk��
     if i(length(i))>=L
         flag=1;
         break;
     else
         if (Rdirection)
             [pk(j),pkp(j)]=max((diffE(i(1):i(length(i)))));
         else
             [pk(j),pkp(j)]=min((diffE(i(1):i(length(i)))));

         end
         pkp(j)=pkp(j)+i(1);        
     end
     
%      if Rdirection
%          i=dermaxx(j)+pk_pre:dermaxx(j)+pk_back;
%      else
%          i=derminx(j)+pk_pre:derminx(j)+pk_back;
%      end
%          if i(length(i))>=L
%             flag=1;
%             break;
%         else
%             pkp1=find(diffE(i(1):i(length(i)-1)).*diffE(i(2):i(length(i)))<=0);
%         end
%         pkp(j)=pkp1(1)+i(1);
%         pk(j)=diffE(pkp(j));
         
     pk(j)=abs(pk(j));
     q_pre=fix(nFreq*200/1000);
     q_b=fix(nFreq*10/1000);
     Qfind=0;
     if Rp(j)>q_pre
         i=Rp(j)-q_b:-1:Rp(j)-q_pre;   %Ѱ��Q����
     else
         i=Rp(j)-q_b:-1:1;
     end
     if i(1)<=0
         break;
     else    
             Qpeakp1=find(diffE(i(2):-1:i(length(i))).*diffE(i(1):-1:i(length(i)-1))<=0);
             if Qpeakp1
                 Qpeakp(j)=i(1)-Qpeakp1(1);
             else
                 Qpeakp(j)=1;
             end
     end  
     Qpeak(j)=FE(Qpeakp(j)); 
     
     t_pre=fix(nFreq*80/1000);
     t_back=fix(nFreq*400/1000);
     i=pkp(j):1:pkp(j)+Swin2 ;    %T����ֵ��λ
     if i(length(i))>=L
         flag=1;
         break;
     else
         if Rdirection
             [Tpeak(j),Tpeakp(j)]=max((FE(i(1):i(length(i)))));  %QRSΪ����ʱ��T���Ϊ�ֲ���Сֵ
         else
             [Tpeak(j),Tpeakp(j)]=min((FE(i(1):i(length(i)))));  %QRSΪ����ʱ��T���Ϊ�ֲ����ֵ
         end
         Tpeakp(j)= Tpeakp(j)+i(1);
     end
       
     tp_pre=fix(nFreq*160/1000);
     if Tpeakp(j)>tp_pre
         i=Tpeakp(j):-1:Tpeakp(j)-tp_pre; %T��ǰ��־���ֵ����
     else
         i=Tpeakp(j):-1:1;
     end
     if i(length(i))<=0
         break;
     else
         [Tdermax(j),Tdermaxp(j)]=max(abs(diffE(i(length(i)):i(1))));
         Tdermaxp(j)=Tdermaxp(j)+i(length(i));
     end
      
%      THTb=Tdermax(j)/2;  
%      t_pre=fix(nFreq*300/1000);
     t_pre=Swin2;
     t_b=fix(nFreq*20/1000);
%      if Tdermaxp(j)>=t_pre
%          i=Tdermaxp(j):-1:Tdermaxp(j)-t_pre;
%      else
%          i=Tdermaxp(j):-1:1;
%      end
     if Tpeakp(j)>t_pre
         i=Tpeakp(j)-t_b:-1:Tpeakp(j)-t_pre; 
     else
         i=Tpeakp(j)-t_b:-1:1;
     end
     if i(length(i))<=0
         break;
     else 
         Tbp1=find(diffE(i(2):-1:i(length(i))).*diffE(i(1):-1:i(length(i)-1))<=0);
         Tbp(j)=i(1)-Tbp1(1);
         Tb(j)=FE(Tbp(j));
%          if Rdirection
%              Tbp1=find((diffE(i(length(i)):i(1)))<=THTb);         
%          else
%              Tbp1=find((diffE(i(length(i)):i(1)))>=(-THTb));
%          end
%          Tbp(j)=Tbp1(length(Tbp1))+i(length(i));
%          Tb(j)=FE(Tbp(j));
     end
    
     tp_back=fix(nFreq*100/1000);
     i=Tpeakp(j):1:Tpeakp(j)+tp_back;  %T����־���ֵ����
     if i(length(i))>=L
         flag=1;
         break;
     else 
         [Tdermin(j),Tderminp(j)]=max(abs(diffE(i(1):i(length(i)))));
         Tderminp(j)=Tderminp(j)+i(1);
     end
       
     THTe=Tdermin(j)/7;        %T���յ���ֵ
     te_back=fix(nFreq*125/1000);
     i=Tderminp(j):1:Tderminp(j)+te_back;
     if i(length(i))>=L
         flag=1;
         break;
     else
         Tep1=find(diffE(i(2):i(length(i))).*diffE(i(1):i(length(i)-1))<=0);
         if Tep1
             Tep(j)=i(1)+Tep1(1);
         else
             Tep(j)=i(length(i));
         end
         Te(j)=FE(Tep(j));
%           if Rdirection
%               Tep1=find(diffE(i(1):i(length(i)))<=(-THTe));
%           else
%               Tep1=find(diffE(i(1):i(length(i)))>=THTe);
%           end
%           if Tep1
%               Tep(j)=Tep1(length(Tep1))+i(1);
% %               Tep(j)=Tep1(1)+i(1);
%               Te(j)=FE(Tep(j));
%           else
%               THTe=THTe/2;
%               break;
%           end
     end 
     
     plot(Rp(j),FE(Rp(j)),'rx')   %��ʾR��λ��
     plot(Qpeakp(j),Qpeak(j),'gx')
     plot(Tpeakp(j),Tpeak(j),'rx')   
     plot(Tbp(j),Tb(j),'mx')   
     plot(Tep(j),Te(j),'mx')
     j=j+1;   
end
% legend('ECG�ź�','R����','Q��','T����','T���','T�յ�');
subplot(2,1,2)
plot(newT,diffE)
hold on
plot(dermaxx,dermay,'r*')
legend('ECG��һ�׵���','�ֲ����ֵ');
title('ECG��һ�׵���');
%  
%%%%%%%%%%%%then find the point of PPG


figure
subplot(3,1,1)
plot(newT,FP)
hold on
ddiffFP=[0,diff(diffFP)]; %second order difference
% ddiffFP=filter(b,a,ddiffFP);

flag=0;
j=1;
while j<=length(Qpeakp) & flag==0
    %%%%%%%%%find the max point B
    if Qpeakp(j)+Swin<L
        i=Qpeakp(j):Qpeakp(j)+Swin;
    else
        flag=1;
        break;
    end
    
    [PPG_B(j),PPG_Bp(j)]=max(FP(i(1):i(length(i))));
    PPG_Bp(j)=PPG_Bp(j)+i(1);
    [pdermay(j),pdermaxx(j)]=max(diffFP(i(1):i(length(i)))); %%һ�׵������ֵ
    pdermaxx(j)=pdermaxx(j)+i(1);  
   
    %%%%%%%%4 methods were used fo find the onset of Pulse wave (A point)
    
    i=Qpeakp(j):PPG_Bp(j); 
    %%%%%%%%%%%%%%%%%%%%%%%%%1. the max second derivative--A
    [Addermax(j),Addermaxp(j)]=max(ddiffFP(i(1):i(length(i))));
    PPG_Ap(j)=Qpeakp(j)+Addermaxp(j);              
    PPG_A(j)=FP(PPG_Ap(j));
    %%%%%%%%%%%%%%%%%%%%%%3.��������������Сֵ��--A3      
    [PPGA3(j),PPGA3p(j)]=min(FP(i(1):i(length(i))));
    PPG_A3p(j)=PPGA3p(j)+i(1);
    PPG_A3(j)=FP(PPG_A3p(j));  
    %%%%%%%%%%%%%%%%%%%%%%%%% 2.use the 10% of the value of B--A2
    THA=(PPG_B(j)-PPG_A3(j))*0.1;       %the value of B is the relative value from the A3 point
    PPGA2=find(FP(i(1):i(length(i)))-PPG_A3(j)<=THA);
    PPG_A2p(j)=i(1)+PPGA2(length(PPGA2));
    PPG_A2(j)=FP(PPG_A2p(j));
    %%%%%%%%%%%%%%%%%%%%%%%4.��intersecting tangents technique ȷ��������A��--A4
    fitN=1;
    InterCorr0=1;
    while InterCorr0>0.999
        fitN=fitN+2;
        InterX=[];
        InterY=[];
        
        for i=1:fitN
            InterX(i)=pdermaxx(j)+i-2;
        end
        InterY=FP(InterX);
        [InterA,S]=polyfit(InterX,InterY,1);
        Yfit=InterA(1)*InterX+InterA(2);
        [InterCorr,InterP]=corrcoef(Yfit,InterY);
        InterCorr0=InterCorr(2);
    end
    
    FitLineX=PPG_A3p(j):pdermaxx(j);
    FitLineY=InterA(1)*FitLineX+InterA(2);
    Above0=find((FitLineY-PPG_A3(j))>0);
    PPG_A4p(j)=Above0(1)+PPG_A3p(j);
    PPG_A4(j)=PPG_A3(j);
%     plot(FitLineX,FitLineY,'c')
%     plot(FitLineX,PPG_A3(j)*ones(length(FitLineX)),'c')
 
    %%%%%%%%%find D point , if there are no D point then give a threshold
    %%%%%%%%%according to 'the derivation wave, first order or second order'
    
 
    %%%%%%%%find C point , first use the second derivation signal-----C---m*
    b_back=fix(nFreq*300/1000);
    i=PPG_Bp(j)+1:PPG_Bp(j)+b_back;
    if i(length(i))>=L
        flag=1;
        break;
    else
        [Addermin2(j),Adderminp2(j)]=min(ddiffFP(i(1):i(length(i))));
        ddmin(j)=i(1)+Adderminp2(j);
% %         ddfirstzero=find(ddiffFP(i(1)+Adderminp2(j):i(length(i)))>=0);
% %         ddfirstzero=find((ddiffFP(i(1)+Adderminp2(j):i(length(i))-1)).*(ddiffFP(i(1)+Adderminp2(j)+1:i(length(i))))<=0);
%         ddfirstzero=find((diffFP(i(1):i(length(i)-1))).*(diffFP(i(2):i(length(i))))<=0);
%         PPG_Cp(j)=i(1)+ddfirstzero(1);
%         PPG_C(j)=FP(PPG_Cp(j));    
    end
%       
    addmin_back=fix(nFreq*100/1000);
    i=ddmin(j):ddmin(j)+addmin_back;
    if i(length(i))>=L          %use the max of the 2nd derivation ----C3
        flag=1;
        break;
    else
        [Addermax2(j),Addermaxp2(j)]=max(ddiffFP(i(1):i(length(i))));  
        PPG_C3p(j)=i(1)+Addermaxp2(j);
        PPG_C3(j)=FP(PPG_C3p(j));    
    end
    
      %%%%%%%%%find C point %%%%%%%%%take use of the T point of ECG---C2
    Te2Tb=Tep(j)-Tbp(j);
    PPG_C2p(j)=PPG_Bp(j)+Te2Tb;
    PPG_C2(j)=FP(PPG_C2p(j));

    plot(PPG_Bp(j),PPG_B(j),'yo')
%     plot(PPG_Cp(j),PPG_C(j),'m*')
    plot(PPG_C2p(j),PPG_C2(j),'m+')
    plot(PPG_C3p(j),PPG_C3(j),'mo')
    plot(PPG_Ap(j),PPG_A(j),'g*') 
    plot(PPG_A2p(j),PPG_A2(j),'g+') 
    plot(PPG_A3p(j),PPG_A3(j),'go')  
    plot(PPG_A4p(j),PPG_A4(j),'gx') 

    j=j+1;
end
legend('PPG','B','C2','C3','A','A2','A3','A4');
title('�������ź�PPG');
subplot(3,1,2)
plot(newT,diffFP)
title('PPG��һ�׵���');

subplot(3,1,3)
plot(newT,ddiffFP)
hold on
plot(PPG_Ap,Addermax,'g*')
hold on
plot(PPG_C3p,Addermax2,'mo')
hold on
plot(PPG_Bp+Adderminp2,Addermin2,'m*')
title('PPG�Ķ��׵���');